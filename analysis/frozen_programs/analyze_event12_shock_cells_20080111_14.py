#!/usr/bin/env python3
"""Search for shock-cell-like radial maxima in the 11--14 Jan 2008 events.

The individual ridge model for every COR2 anchor is selected exactly as in the
event phase-jitter test.  Along that selected travel-time path, minimally
smoothed robust BFF profiles are constructed independently in tB and pB.
Peaks must be present in both products within 0.10 R_sun.  The complete
product/model selection and cell count are repeated at unrelated times.

This is an axial compression-cell diagnostic, not a proof of two-dimensional
shock diamonds or of MHD critical points.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks

import analyze_pds_nonlinear_transport_20080111_14 as core
import analyze_pds_event_phase_jitter_20080111_14 as event_test


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "results" / "pds_20080111_14_nonlinear_level0"
PHASE_SOURCE = SOURCE / "event_phase_jitter"
OUT = SOURCE / "event12_shock_cells"
OUT.mkdir(parents=True, exist_ok=True)

CONFIGS = {
    "moderate": {"sigma_samples": 1.0, "prominence": 0.50},
    "strict": {"sigma_samples": 1.0, "prominence": 0.75},
    "coarse": {"sigma_samples": 2.0, "prominence": 0.50},
}
MIN_PEAK_DISTANCE_SAMPLES = 3
MATCH_TOLERANCE_RSUN = 0.10


def robust_standardize(raw):
    center = np.nanmedian(raw, axis=0)
    mad = 1.4826 * np.nanmedian(np.abs(raw - center[None]), axis=0)
    std = np.nanstd(raw, axis=0)
    scale = np.where(mad > 1e-8, mad, np.where(std > 1e-8, std, 1.0))
    return np.clip((raw - center[None]) / scale[None], -6.0, 6.0)


def sample_profile(zmap, grid, anchor, tau):
    return np.array([
        np.interp(anchor + tau[j], grid, zmap[:, j], left=np.nan, right=np.nan)
        for j in range(len(tau))
    ])


def detect(profile, sigma_samples, prominence):
    smooth = (
        gaussian_filter1d(np.asarray(profile, float), sigma_samples, mode="nearest")
        if sigma_samples > 0
        else np.asarray(profile, dtype=float).copy()
    )
    peaks, properties = find_peaks(
        smooth,
        prominence=prominence,
        distance=MIN_PEAK_DISTANCE_SAMPLES,
    )
    return smooth, peaks, properties


def match_peaks(radii, primary, secondary):
    candidates = []
    for i in primary:
        for j in secondary:
            distance = abs(float(radii[i] - radii[j]))
            if distance <= MATCH_TOLERANCE_RSUN + 1e-9:
                candidates.append((distance, int(i), int(j)))
    used_primary, used_secondary, matches = set(), set(), []
    for distance, i, j in sorted(candidates):
        if i in used_primary or j in used_secondary:
            continue
        used_primary.add(i)
        used_secondary.add(j)
        matches.append((i, j, distance))
    return sorted(matches)


def metrics_for_config(radii, profiles, primary_product, config):
    detected = {}
    for product in ("total_b", "pb"):
        smooth, peaks, props = detect(
            profiles[product], config["sigma_samples"], config["prominence"]
        )
        detected[product] = {"smooth": smooth, "peaks": peaks, "props": props}
    secondary_product = "pb" if primary_product == "total_b" else "total_b"
    matches = match_peaks(
        radii,
        detected[primary_product]["peaks"],
        detected[secondary_product]["peaks"],
    )
    primary_indices = np.array([item[0] for item in matches], dtype=int)
    secondary_indices = np.array([item[1] for item in matches], dtype=int)
    primary_radii = radii[primary_indices] if len(primary_indices) else np.array([])
    matched_centers = (
        0.5 * (radii[primary_indices] + radii[secondary_indices])
        if len(primary_indices) else np.array([])
    )
    spacing = np.diff(primary_radii)
    spacing_cv = (
        float(np.std(spacing) / np.mean(spacing))
        if len(spacing) >= 2 and np.mean(spacing) > 0 else np.nan
    )
    regularity = 1.0 / (1.0 + spacing_cv) if np.isfinite(spacing_cv) else 0.0
    count = int(len(matches))
    cell_stat = float(count + regularity) if count >= 3 else float(count)
    return {
        "primary_product": primary_product,
        "secondary_product": secondary_product,
        "primary_peak_count": int(len(detected[primary_product]["peaks"])),
        "secondary_peak_count": int(len(detected[secondary_product]["peaks"])),
        "matched_peak_count": count,
        "matched_primary_radii_rsun": primary_radii.tolist(),
        "matched_center_radii_rsun": matched_centers.tolist(),
        "spacing_rsun": spacing.tolist(),
        "spacing_cv": None if not np.isfinite(spacing_cv) else spacing_cv,
        "cell_stat": cell_stat,
        "detected": detected,
        "matches": matches,
    }


def analyze_anchor(anchor, grid, radii, zmaps, raw_zmaps, models, tau_matrix):
    best_score, product, model_index, model = event_test.best_scanned_model(
        zmaps, grid, float(anchor), models, tau_matrix
    )
    profiles = {
        name: sample_profile(zmap, grid, float(anchor), model.tau)
        for name, zmap in raw_zmaps.items()
    }
    configs = {
        name: metrics_for_config(radii, profiles, product, config)
        for name, config in CONFIGS.items()
    }
    return {
        "best_score": float(best_score),
        "best_product": product,
        "best_model": event_test.model_record(model),
        "profiles": profiles,
        "configs": configs,
    }


def empirical_p(observed, null):
    null = np.asarray(null, float)
    return float((np.sum(null >= observed) + 1) / (len(null) + 1))


def stripped_config(row):
    return {key: value for key, value in row.items() if key not in ("detected", "matches")}


def make_figure(event, radii, z_published, null_counts, null_stats):
    strict = event["configs"]["strict"]
    moderate = event["configs"]["moderate"]
    coarse = event["configs"]["coarse"]
    fig, axes = plt.subplots(3, 1, figsize=(12.5, 10.5), constrained_layout=True)

    colors = {"total_b": "#a23b72", "pb": "#4c78a8"}
    for product in ("total_b", "pb"):
        d = strict["detected"][product]
        axes[0].plot(radii, d["smooth"], lw=1.8, color=colors[product], label=product)
        axes[0].scatter(radii[d["peaks"]], d["smooth"][d["peaks"]], s=38,
                        color=colors[product], edgecolor="white", zorder=3)
    for i, j, _ in strict["matches"]:
        x1, x2 = radii[i], radii[j]
        y1 = strict["detected"][strict["primary_product"]]["smooth"][i]
        y2 = strict["detected"][strict["secondary_product"]]["smooth"][j]
        axes[0].plot([x1, x2], [y1, y2], color="#2a9d8f", lw=2.5, alpha=0.8)
    axes[0].axvline(event["best_model"]["change_radius_rsun"], color="#e09f3e", ls="--",
                    label="best-model break")
    axes[0].set(title=(
        f"Event #{event['event_number']}: strict minimally smoothed profiles; "
        f"{strict['matched_peak_count']} matched tB/pB maxima"),
        ylabel="Robust BFF residual z")
    axes[0].legend(ncol=3)
    axes[0].grid(alpha=0.2)

    axes[1].plot(radii, z_published["total_b"], color=colors["total_b"], lw=1.8, label="tB")
    axes[1].plot(radii, z_published["pb"], color=colors["pb"], lw=1.5, label="pB")
    axes[1].axvline(event["best_model"]["change_radius_rsun"], color="#e09f3e", ls="--")
    axes[1].set(title="Published-comparison 13-sample radial smoothing",
                ylabel="Radially smoothed BFF z")
    axes[1].legend()
    axes[1].grid(alpha=0.2)

    bins = np.arange(-0.5, max(6.5, np.max(null_counts) + 1.5), 1.0)
    axes[2].hist(null_counts, bins=bins, color="0.72", edgecolor="white",
                 label="unrelated-time null")
    axes[2].axvline(strict["matched_peak_count"], color="#a23b72", lw=2.5,
                    label=f"event #{event['event_number']} count={strict['matched_peak_count']}")
    axes[2].set(title=(
        f"Strict matched-cell count null: p={event['strict_count_p_raw']:.3f}; "
        f"count+regularity p={event['strict_cell_stat_p_raw']:.3f}"),
        xlabel="Matched tB/pB maxima", ylabel="Null anchors")
    axes[2].legend()
    axes[2].grid(axis="y", alpha=0.2)
    for ax in axes[:2]:
        ax.set_xlabel(r"Radius ($R_\odot$)")
    stamp = event["utc"].replace("-", "").replace(":", "").replace("T", "_")[:13]
    path = OUT / f"pds_{stamp}_event{event['event_number']:02d}_shock_cell_diagnostic.png"
    fig.savefig(path, dpi=230)
    plt.close(fig)
    return path


def write_report(
    rows, event12, event9, null_count, published_counts12, published_counts9,
    event9_scan_excursion_probability,
):
    strict = event12["configs"]["strict"]
    moderate = event12["configs"]["moderate"]
    coarse = event12["configs"]["coarse"]
    lines = [
        "# 11--14 January 2008 shock-cell / local-maxima test",
        "",
        "## Verdict",
        "",
        f"At the predeclared strict minimal-smoothing scale, event #12 contains "
        f"{strict['matched_peak_count']} tB/pB-matched radial maxima at "
        + ", ".join(f"{value:.3f} R_sun" for value in strict["matched_center_radii_rsun"]) + ".",
        f"However, the complete unrelated-time null gives raw p={event12['strict_count_p_raw']:.4f} "
        f"for the matched count and p={event12['strict_cell_stat_p_raw']:.4f} after also rewarding "
        "regular spacing.  Therefore the three-maxima pattern is not yet an unusual shock-cell detection.",
        "",
        f"The count is scale dependent: moderate settings give {moderate['matched_peak_count']} matched "
        f"maxima, strict settings give {strict['matched_peak_count']}, and coarse smoothing gives "
        f"{coarse['matched_peak_count']}.  With the 13-sample radial smoothing used in the published "
        f"comparison map, the tB/pB profiles contain {published_counts12['total_b']}/"
        f"{published_counts12['pb']} prominent maxima separately, not a three-cell common train.",
        "",
        "## Strongest local morphology candidate: event #9",
        "",
        f"Event #9 (2008-01-13 15:37 UT) has "
        f"{event9['configs']['strict']['matched_peak_count']} strict matched maxima at "
        + ", ".join(
            f"{value:.3f} R_sun"
            for value in event9["configs"]["strict"]["matched_center_radii_rsun"]
        ) + ".",
        f"Its spacing CV is {event9['configs']['strict']['spacing_cv']:.3f}; the raw local "
        f"count+regularity p value is {event9['strict_cell_stat_p_raw']:.4f}.  It is more "
        f"scale-stable than event #12: the moderate/strict/coarse matched counts are "
        f"{event9['configs']['moderate']['matched_peak_count']}/"
        f"{event9['configs']['strict']['matched_peak_count']}/"
        f"{event9['configs']['coarse']['matched_peak_count']}, while the published-smoothed "
        f"tB/pB profiles contain {published_counts9['total_b']}/{published_counts9['pb']} "
        f"prominent maxima separately and {published_counts9['matched']} matched maxima.",
        f"This is still not a discovery after choosing the best morphology among 12 events.  "
        f"Using the unrelated-time null directly, the probability that a 12-anchor scan contains "
        f"at least one cell statistic as large as event #9 is about "
        f"{event9_scan_excursion_probability:.3f}.  This is a selection-aware max-excursion check, "
        "not a BH correction.",
        "The independently tested propagating ridge for event #9 is also not significant "
        "(raw p=0.486), and its phase kink is not identifiable (raw p=0.819).  Thus event #9 is "
        "a useful morphology candidate, not yet a causal shock-cell detection.",
        "",
        "## Habbal critical points are not image maxima",
        "",
        "In the multiple-transonic-solution models of Habbal and collaborators, three critical points "
        "are singular/regularity points of the steady wind equation generated by localized momentum "
        "addition or rapid flow-tube divergence.  The allowed global solution may pass continuously "
        "through one critical point or connect critical solution branches through a standing shock.  "
        "Three brightness maxima therefore do not represent three critical points.",
        "",
        "A shock-diamond interpretation would require a two-dimensional sequence of alternating oblique "
        "fronts/X nodes, approximately repeatable cell spacing, tB/pB agreement, and persistence under "
        "reasonable filtering.  The present one-dimensional axial test finds a suggestive three-maximum "
        "configuration at one processing scale, but it fails the null and scale-robustness requirements.",
        "",
        "## Event-by-event strict results",
        "",
        "| # | UTC | Best product/model | Matched maxima | Spacing CV | Count p raw | Cell-stat p raw |",
        "| ---: | --- | --- | ---: | ---: | ---: | ---: |",
    ]
    for row in rows:
        model = row["best_model"]
        if model["change_radius_rsun"] is None:
            model_text = f"{row['best_product']} {model['family']} {model['v_inner_km_s']:.0f}"
        else:
            model_text = (
                f"{row['best_product']} {model['family']} "
                f"{model['v_inner_km_s']:.0f}->{model['v_outer_km_s']:.0f} @"
                f"{model['change_radius_rsun']:.2f}"
            )
        cv = row["configs"]["strict"]["spacing_cv"]
        lines.append(
            f"| {row['event_number']} | {row['utc'].replace('T', ' ')[:16]} | {model_text} | "
            f"{row['configs']['strict']['matched_peak_count']} | "
            f"{'--' if cv is None else f'{cv:.3f}'} | {row['strict_count_p_raw']:.3f} | "
            f"{row['strict_cell_stat_p_raw']:.3f} |"
        )
    lines += [
        "",
        "## Method",
        "",
        "- The event-specific product and constant/acceleration/deceleration model are selected by the "
        "previous full scan.",
        "- Unsmeared temporal BFF residuals are robust-standardized at every height, then smoothed by only "
        "one 0.025 R_sun radial sample for the strict test.",
        "- Peaks require prominence >=0.75 and separation >=0.075 R_sun; tB and pB peaks must agree "
        "within 0.10 R_sun.",
        f"- The complete selection and cell measurement are repeated at {null_count} unrelated times.",
        "- All p values are raw; no BH adjustment.",
        "",
        "## Limitation",
        "",
        "The COR1 products are direct bias/exposure-normalized Level-0 polarization fits rather than a "
        "complete SECCHI_PREP/background calibration.  The test cannot provide a density compression "
        "ratio or magnetosonic Mach number.",
    ]
    path = OUT / "PDS_20080111_14_shock_cells_no_BH.md"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def main():
    maps_path = SOURCE / "pds_20080111_14_level0_nonlinear_transport_maps.npz"
    with np.load(maps_path) as data:
        grid = np.asarray(data["minutes"], float)
        radii = np.asarray(data["radii"], float)
        zmaps = {
            "total_b": np.asarray(data["z_total_b"], float),
            "pb": np.asarray(data["z_pb"], float),
        }
        raw_maps = {
            "total_b": np.asarray(data["raw_total_b"], float),
            "pb": np.asarray(data["raw_pb"], float),
        }
        event_times = np.asarray(data["event_times"], float)
    raw_zmaps = {name: robust_standardize(raw) for name, raw in raw_maps.items()}

    phase_results = json.loads(
        (PHASE_SOURCE / "pds_20080111_14_event_phase_jitter_results.json").read_text()
    )
    core.RADII_WORK = radii
    families = core.model_families(radii)
    models = []
    for family in ("constant", "acceleration", "deceleration"):
        models.extend(families[family])
    tau_matrix = np.stack([model.tau for model in models])

    reference = event_test.model_from_dict(phase_results["reference_model"], radii)
    margin = abs(event_test.PHASE_GRID[0]) + abs(event_test.LOCAL_GRID[0])
    lo = grid[0] - np.min(reference.tau) + margin
    hi = grid[-1] - margin
    field_anchors = np.arange(math.ceil(lo / 30.0) * 30.0,
                              math.floor(hi / 30.0) * 30.0 + 0.1, 30.0)

    null_rows = []
    for number, anchor in enumerate(field_anchors, 1):
        null_rows.append(analyze_anchor(
            float(anchor), grid, radii, zmaps, raw_zmaps, models, tau_matrix
        ))
        if number % 50 == 0:
            print(f"Cell null: {number}/{len(field_anchors)}", flush=True)
    null_counts = np.array([
        row["configs"]["strict"]["matched_peak_count"] for row in null_rows
    ])
    null_stats = np.array([
        row["configs"]["strict"]["cell_stat"] for row in null_rows
    ])

    rows = []
    full_events = []
    for i, anchor in enumerate(event_times):
        event = analyze_anchor(float(anchor), grid, radii, zmaps, raw_zmaps, models, tau_matrix)
        strict = event["configs"]["strict"]
        event["event_number"] = i + 1
        event["utc"] = core.STRONG_EVENTS[i][0]
        event["cor2_z"] = float(core.STRONG_EVENTS[i][1])
        event["strict_count_p_raw"] = empirical_p(strict["matched_peak_count"], null_counts)
        event["strict_cell_stat_p_raw"] = empirical_p(strict["cell_stat"], null_stats)
        full_events.append(event)
        rows.append({
            "event_number": event["event_number"],
            "utc": event["utc"],
            "cor2_z": event["cor2_z"],
            "best_score": event["best_score"],
            "best_product": event["best_product"],
            "best_model": event["best_model"],
            "strict_count_p_raw": event["strict_count_p_raw"],
            "strict_cell_stat_p_raw": event["strict_cell_stat_p_raw"],
            "configs": {name: stripped_config(value) for name, value in event["configs"].items()},
        })

    event9 = full_events[8]
    event12 = full_events[11]

    def published_diagnostic(event, event_index):
        selected_model = event_test.model_from_dict(event["best_model"], radii)
        profiles = {
            name: sample_profile(zmap, grid, float(event_times[event_index]), selected_model.tau)
            for name, zmap in zmaps.items()
        }
        counts = {}
        peak_indices = {}
        for name, profile in profiles.items():
            _, peaks, _ = detect(profile, 0.0, 0.30)
            counts[name] = int(len(peaks))
            peak_indices[name] = peaks
        counts["matched"] = int(len(match_peaks(
            radii, peak_indices["total_b"], peak_indices["pb"]
        )))
        return profiles, counts

    published_profiles9, published_counts9 = published_diagnostic(event9, 8)
    published_profiles12, published_counts12 = published_diagnostic(event12, 11)

    event9_exceedances = int(np.sum(null_stats >= event9["configs"]["strict"]["cell_stat"]))
    scan_size = len(full_events)
    event9_scan_excursion_probability = float(
        1.0
        - math.comb(len(null_stats) - event9_exceedances, scan_size)
        / math.comb(len(null_stats), scan_size)
    )

    figure_path9 = make_figure(event9, radii, published_profiles9, null_counts, null_stats)
    figure_path12 = make_figure(event12, radii, published_profiles12, null_counts, null_stats)
    report_path = write_report(
        rows, event12, event9, len(field_anchors), published_counts12, published_counts9,
        event9_scan_excursion_probability,
    )

    csv_path = OUT / "pds_20080111_14_shock_cell_event_table.csv"
    flat_rows = []
    for row in rows:
        strict = row["configs"]["strict"]
        flat_rows.append({
            "event_number": row["event_number"],
            "utc": row["utc"],
            "cor2_z": row["cor2_z"],
            "best_product": row["best_product"],
            "best_family": row["best_model"]["family"],
            "best_v_inner_km_s": row["best_model"]["v_inner_km_s"],
            "best_v_outer_km_s": row["best_model"]["v_outer_km_s"],
            "best_change_radius_rsun": row["best_model"]["change_radius_rsun"],
            "strict_matched_peak_count": strict["matched_peak_count"],
            "strict_matched_center_radii_rsun": "|".join(
                f"{x:.3f}" for x in strict["matched_center_radii_rsun"]
            ),
            "strict_spacing_cv": strict["spacing_cv"],
            "strict_count_p_raw": row["strict_count_p_raw"],
            "strict_cell_stat_p_raw": row["strict_cell_stat_p_raw"],
        })
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(flat_rows[0]))
        writer.writeheader()
        writer.writerows(flat_rows)

    output = {
        "event": "2008-01-11--14",
        "target_event_number": 12,
        "target_event_utc": core.STRONG_EVENTS[11][0],
        "criteria": {
            **CONFIGS,
            "minimum_peak_distance_samples": MIN_PEAK_DISTANCE_SAMPLES,
            "minimum_peak_distance_rsun": MIN_PEAK_DISTANCE_SAMPLES * float(np.median(np.diff(radii))),
            "tb_pb_match_tolerance_rsun": MATCH_TOLERANCE_RSUN,
        },
        "field_null_count": int(len(field_anchors)),
        "strongest_local_morphology_event_number": 9,
        "event9": rows[8],
        "event9_scan_excursion_probability_12_anchors": event9_scan_excursion_probability,
        "event12": rows[11],
        "events": rows,
        "null_strict_matched_counts": null_counts.tolist(),
        "null_strict_cell_stats": null_stats.tolist(),
        "published_smoothing_peak_counts_event9": published_counts9,
        "published_smoothing_peak_counts_event12": published_counts12,
        "no_BH": True,
        "files": {
            "report": report_path.name,
            "figure_event9": figure_path9.name,
            "figure_event12": figure_path12.name,
            "table": csv_path.name,
        },
    }
    json_path = OUT / "pds_20080113_2137_event12_shock_cells_results.json"
    json_path.write_text(json.dumps(output, indent=2), encoding="utf-8")
    print(json.dumps({
        "event12_configs": rows[11]["configs"],
        "event12_count_p_raw": rows[11]["strict_count_p_raw"],
        "event12_cell_stat_p_raw": rows[11]["strict_cell_stat_p_raw"],
        "event9_configs": rows[8]["configs"],
        "event9_cell_stat_p_raw": rows[8]["strict_cell_stat_p_raw"],
        "event9_scan_excursion_probability_12_anchors": event9_scan_excursion_probability,
        "all_event_strict_counts": [row["configs"]["strict"]["matched_peak_count"] for row in rows],
        "all_event_count_p_raw": [row["strict_count_p_raw"] for row in rows],
        "published_counts_event9": published_counts9,
        "published_counts_event12": published_counts12,
    }, indent=2))


if __name__ == "__main__":
    main()
