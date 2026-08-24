#!/usr/bin/env python3
"""Frozen alternative-representation robustness test for events #9 and #12."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import analyze_event09_2d_xfront_20080113 as x2d
import analyze_pds_event_phase_jitter_20080111_14 as phase
import analyze_pds_nonlinear_transport_20080111_14 as core
import analyze_two_candidate_2d_xfront_20080113 as pair_test


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "results" / "pds_20080111_14_nonlinear_level0"
OUT = SOURCE / "two_event_representation_robustness"
OUT.mkdir(parents=True, exist_ok=True)

PAIR_NULL_FILE = (
    SOURCE / "two_candidate_2d_xfront"
    / "pds_20080113_event09_event12_fixed_360min_pair_null.csv"
)

CANDIDATES = {
    "event9": {
        "event_number": 9,
        "event_index": 8,
        "utc": core.STRONG_EVENTS[8][0],
        "model": {
            "family": "deceleration",
            "v_inner_km_s": 100.0,
            "v_outer_km_s": 30.0,
            "change_radius_rsun": 2.20,
        },
        "nodes": np.array([1.775, 2.125, 2.500, 2.825], float),
        "required_positive": 3,
    },
    "event12": {
        "event_number": 12,
        "event_index": 11,
        "utc": core.STRONG_EVENTS[11][0],
        "model": {
            "family": "deceleration",
            "v_inner_km_s": 100.0,
            "v_outer_km_s": 25.0,
            "change_radius_rsun": 2.90,
        },
        "nodes": np.array([1.950, 2.350, 2.700], float),
        "required_positive": 3,
    },
}

REPRESENTATIONS = ("base60", "base120", "nrgf60")
PAIR_SEPARATION_MIN = 360.0


def robust_temporal_standardize(values):
    values = np.asarray(values, float)
    center = np.nanmedian(values, axis=0)
    mad = 1.4826 * np.nanmedian(np.abs(values - center[None]), axis=0)
    std = np.nanstd(values, axis=0)
    scale = np.where(mad > 1e-8, mad, np.where(std > 1e-8, std, 1.0))
    return np.clip((values - center[None]) / scale[None], -6.0, 6.0)


def lag_difference(values, frames):
    output = np.full_like(values, np.nan, dtype=float)
    output[frames:] = values[frames:] - values[:-frames]
    return output


def nrgf(values):
    mean = np.nanmean(values, axis=2, keepdims=True)
    std = np.nanstd(values, axis=2, keepdims=True)
    std = np.where(std > 1e-8, std, 1.0)
    return (values - mean) / std


def regularize_cube():
    with np.load(x2d.CUBE_FILE) as data:
        cubes = {name: np.asarray(data[name], float) for name in ("total_b", "pb")}
        minutes = np.asarray(data["minutes"], float)
        radii_full = np.asarray(data["radii"], float)
        pa = np.asarray(data["pa"], float)
        quality = np.asarray(data["fit_quality"], float)
    with np.load(x2d.MAPS_FILE) as data:
        radii = np.asarray(data["radii"], float)
        path_full = np.asarray(data["traced_path_pa"], float)
        event_times = np.asarray(data["event_times"], float)

    center = np.nanmedian(quality)
    mad = 1.4826 * np.nanmedian(np.abs(quality - center))
    good = np.isfinite(quality) & (quality <= center + 6 * max(mad, 1e-8))
    grid = np.arange(
        math.ceil(np.min(minutes[good]) / core.CADENCE_MIN) * core.CADENCE_MIN,
        math.floor(np.max(minutes[good]) / core.CADENCE_MIN) * core.CADENCE_MIN + 0.1,
        core.CADENCE_MIN,
    )
    regular = {}
    for product, cube in cubes.items():
        flat = cube.reshape(len(cube), -1)
        output = np.empty((len(grid), flat.shape[1]), float)
        for j in range(flat.shape[1]):
            valid = good & np.isfinite(flat[:, j])
            output[:, j] = (
                np.interp(grid, minutes[valid], flat[valid, j]) if np.sum(valid) >= 2 else np.nan
            )
        regular[product] = output.reshape((len(grid),) + cube.shape[1:])
    return grid, radii_full, radii, pa, path_full, event_times, regular


def build_representations(regular, radii_full, radii, pa, path_full):
    cadence = core.CADENCE_MIN
    lag60 = int(round(60.0 / cadence))
    lag120 = int(round(120.0 / cadence))
    output = {name: {} for name in REPRESENTATIONS}
    for product, cube in regular.items():
        derived = {
            "base60": lag_difference(cube, lag60),
            "base120": lag_difference(cube, lag120),
            "nrgf60": lag_difference(nrgf(cube), lag60),
        }
        for name, residual in derived.items():
            zcube = robust_temporal_standardize(residual)
            output[name][product] = x2d.path_offset_cube(
                zcube, radii_full, pa, path_full, radii
            )
            del zcube
    return output


def event_maps(offset_products, grid, anchor, model):
    return {
        product: x2d.aligned_slice(cube, grid, anchor, model.tau)
        for product, cube in offset_products.items()
    }


def score_candidate(offset_products, grid, anchor, model, radii, nodes):
    maps = event_maps(offset_products, grid, anchor, model)
    score = pair_test.score_top_n(maps, radii, nodes, len(nodes))
    return maps, score


def empirical_p(exceedances, total):
    return float((int(exceedances) + 1) / (int(total) + 1))


def coherent_map(maps):
    same = maps["total_b"] * maps["pb"] > 0
    result = np.zeros_like(maps["total_b"])
    result[same] = np.sign(maps["total_b"][same]) * np.sqrt(
        np.abs(maps["total_b"][same] * maps["pb"][same])
    )
    return result


def make_maps_figure(event_results, radii):
    fig, axes = plt.subplots(3, 2, figsize=(13.5, 14), constrained_layout=True)
    labels = {
        "base60": "Base difference 60 min",
        "base120": "Base difference 120 min",
        "nrgf60": "NRGF difference 60 min",
    }
    for row, representation in enumerate(REPRESENTATIONS):
        for col, candidate_name in enumerate(("event9", "event12")):
            candidate = CANDIDATES[candidate_name]
            result = event_results[candidate_name][representation]
            image = axes[row, col].pcolormesh(
                x2d.PA_OFFSETS, radii, coherent_map(result["maps"]),
                cmap="RdBu_r", vmin=-2.5, vmax=2.5, shading="auto",
            )
            axes[row, col].axvline(0, color="0.2", lw=0.9, ls=":")
            for node, node_fit in zip(candidate["nodes"], result["score"]["nodes"]):
                axes[row, col].axhline(node, color="black", lw=0.7, ls=":")
                q = np.linspace(-8, 8, 81)
                vertex = node_fit["vertex_radius_rsun"]
                slope = node_fit["absolute_slope_rsun_per_degree"]
                axes[row, col].plot(q, vertex + slope * q, color="#f6bd60", lw=0.9)
                axes[row, col].plot(q, vertex - slope * q, color="#84a59d", lw=0.9)
            axes[row, col].set(
                title=(f"{labels[representation]} | event #{candidate['event_number']} | "
                       f"median X={result['score']['median_coherent_x_response']:.3f}"),
                xlabel="PA offset (deg)", ylabel=r"Radius ($R_\odot$)",
            )
            fig.colorbar(image, ax=axes[row, col], label="Same-sign tB/pB z")
    path = OUT / "pds_20080113_two_event_representation_maps.png"
    fig.savefig(path, dpi=230)
    plt.close(fig)
    return path


def make_null_figure(pair_rows, event_results, summary):
    fig, axes = plt.subplots(2, 2, figsize=(13, 10), constrained_layout=True)
    colors = {"base60": "#4c78a8", "base120": "#e45756", "nrgf60": "#72b7b2"}
    for ax, representation in zip(axes.flat[:3], REPRESENTATIONS):
        early = np.array([row[f"{representation}_early_score"] for row in pair_rows], float)
        late = np.array([row[f"{representation}_late_score"] for row in pair_rows], float)
        valid = np.isfinite(early) & np.isfinite(late)
        threshold_early = event_results["event9"][representation]["score"][
            "median_coherent_x_response"
        ]
        threshold_late = event_results["event12"][representation]["score"][
            "median_coherent_x_response"
        ]
        ax.scatter(early[valid], late[valid], s=18, alpha=0.55, color=colors[representation])
        ax.axvline(threshold_early, color="black", ls="--")
        ax.axhline(threshold_late, color="black", ls="--")
        ax.scatter([threshold_early], [threshold_late], marker="*", s=150,
                   color="#9b2c65", edgecolor="white", zorder=4, label="real pair")
        ax.set(
            title=f"{representation}: pair p={summary['representations'][representation]['p_raw']:.3f}",
            xlabel="Early (#9-like) median X", ylabel="Late (#12-like) median X",
        )
        ax.legend()
        ax.grid(alpha=0.2)

    counts = np.array([row["representation_exceedance_count"] for row in pair_rows], int)
    axes[1, 1].hist(counts, bins=np.arange(-0.5, 4.5, 1), color="0.72", edgecolor="white")
    axes[1, 1].axvline(2, color="#9b2c65", lw=2.5, label="primary >=2 representations")
    axes[1, 1].set(
        title=(f"Multi-representation pair null: "
               f"p={summary['primary_pair_p_raw']:.4f}"),
        xlabel="Representations exceeding both real-event scores", ylabel="Shifted pairs",
    )
    axes[1, 1].legend()
    axes[1, 1].grid(axis="y", alpha=0.2)
    path = OUT / "pds_20080113_two_event_representation_nulls.png"
    fig.savefig(path, dpi=230)
    plt.close(fig)
    return path


def write_report(summary):
    lines = [
        "# Two-event alternative-representation robustness test",
        "",
        "## Verdict",
        "",
        summary["verdict"],
        "",
        "## Real-event scores and fixed-pair null (raw, no BH)",
        "",
        "| Representation | #9 median X | #9 positive | #12 median X | #12 positive | "
        "Pair exceedances | Pair p |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for representation in REPRESENTATIONS:
        row = summary["representations"][representation]
        lines.append(
            f"| {representation} | {row['event9_median_x']:.4f} | "
            f"{row['event9_positive_nodes']}/4 | {row['event12_median_x']:.4f} | "
            f"{row['event12_positive_nodes']}/3 | {row['pair_exceedances']}/"
            f"{row['valid_pairs']} | {row['p_raw']:.4f} |"
        )
    lines += [
        "",
        "## Primary >=2-of-3 result",
        "",
        f"- Valid shifted pairs: {summary['valid_pair_count']}.",
        f"- Shifted pairs exceeding the real pair in at least two representations: "
        f"{summary['primary_pair_exceedances']}.",
        f"- Raw primary pair p={summary['primary_pair_p_raw']:.4f}.",
        f"- Event #9 meets its >=75% positive-node condition in "
        f"{summary['event9_positive_representation_count']}/3 representations.",
        f"- Event #12 meets its >=75% positive-node condition in "
        f"{summary['event12_positive_representation_count']}/3 representations.",
        "",
        "## Timing and geometry",
        "",
        "No global phase shift, height-dependent timing correction, PA shift, node reselection, "
        "or separation scan was fitted.  The real #9 -> #12 separation remains exactly 360 min.",
        "",
        "## Meaning and limitation",
        "",
        "This test asks whether the ordered pair survives alternative image transformations and "
        "therefore addresses sensitivity to the original BFF filter.  The three representations "
        "are correlated transformations of the same Level-0-derived measurements; they are not "
        "three independent data sets.  SECCHI_PREP Level-1/background validation and a physical "
        "density-jump/Mach-number test remain necessary for a shock claim.",
    ]
    path = OUT / "PDS_20080113_two_event_representation_robustness_no_BH.md"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def main():
    grid, radii_full, radii, pa, path_full, event_times, regular = regularize_cube()
    representations = build_representations(regular, radii_full, radii, pa, path_full)

    models = {
        name: phase.model_from_dict(candidate["model"], radii)
        for name, candidate in CANDIDATES.items()
    }
    event_results = {name: {} for name in CANDIDATES}
    for candidate_name, candidate in CANDIDATES.items():
        anchor = float(event_times[candidate["event_index"]])
        for representation in REPRESENTATIONS:
            maps, score = score_candidate(
                representations[representation], grid, anchor, models[candidate_name],
                radii, candidate["nodes"],
            )
            event_results[candidate_name][representation] = {"maps": maps, "score": score}

    with PAIR_NULL_FILE.open(newline="", encoding="utf-8") as handle:
        early_anchors = np.array(
            [float(row["early_anchor_min"]) for row in csv.DictReader(handle)], float
        )

    pair_rows = []
    for number, early_anchor in enumerate(early_anchors, 1):
        late_anchor = float(early_anchor + PAIR_SEPARATION_MIN)
        row = {"early_anchor_min": float(early_anchor), "late_anchor_min": late_anchor}
        fully_valid = True
        exceedance_count = 0
        for representation in REPRESENTATIONS:
            _, early_score = score_candidate(
                representations[representation], grid, float(early_anchor), models["event9"],
                radii, CANDIDATES["event9"]["nodes"],
            )
            _, late_score = score_candidate(
                representations[representation], grid, late_anchor, models["event12"],
                radii, CANDIDATES["event12"]["nodes"],
            )
            early_value = early_score["median_coherent_x_response"]
            late_value = late_score["median_coherent_x_response"]
            early_threshold = event_results["event9"][representation]["score"][
                "median_coherent_x_response"
            ]
            late_threshold = event_results["event12"][representation]["score"][
                "median_coherent_x_response"
            ]
            valid = bool(np.isfinite(early_value) and np.isfinite(late_value))
            passed = bool(
                valid and early_value >= early_threshold and late_value >= late_threshold
            )
            row[f"{representation}_early_score"] = float(early_value)
            row[f"{representation}_late_score"] = float(late_value)
            row[f"{representation}_pair_pass"] = passed
            fully_valid &= valid
            exceedance_count += int(passed)
        row["fully_valid"] = fully_valid
        row["representation_exceedance_count"] = exceedance_count
        row["primary_pair_pass"] = bool(fully_valid and exceedance_count >= 2)
        pair_rows.append(row)
        if number % 25 == 0:
            print(f"Representation pair null: {number}/{len(early_anchors)}", flush=True)

    valid_rows = [row for row in pair_rows if row["fully_valid"]]
    representation_summary = {}
    for representation in REPRESENTATIONS:
        exceedances = int(sum(row[f"{representation}_pair_pass"] for row in valid_rows))
        e9score = event_results["event9"][representation]["score"]
        e12score = event_results["event12"][representation]["score"]
        representation_summary[representation] = {
            "event9_median_x": float(e9score["median_coherent_x_response"]),
            "event9_positive_nodes": int(e9score["positive_node_count"]),
            "event12_median_x": float(e12score["median_coherent_x_response"]),
            "event12_positive_nodes": int(e12score["positive_node_count"]),
            "pair_exceedances": exceedances,
            "valid_pairs": int(len(valid_rows)),
            "p_raw": empirical_p(exceedances, len(valid_rows)),
        }

    primary_exceedances = int(sum(row["primary_pair_pass"] for row in valid_rows))
    primary_p = empirical_p(primary_exceedances, len(valid_rows))
    event9_positive_reps = int(sum(
        representation_summary[name]["event9_positive_nodes"]
        >= CANDIDATES["event9"]["required_positive"]
        for name in REPRESENTATIONS
    ))
    event12_positive_reps = int(sum(
        representation_summary[name]["event12_positive_nodes"]
        >= CANDIDATES["event12"]["required_positive"]
        for name in REPRESENTATIONS
    ))
    passed = bool(
        primary_p < 0.05 and event9_positive_reps >= 2 and event12_positive_reps >= 2
    )
    if passed:
        verdict = (
            "The ordered #9 -> #12 pair passes the frozen alternative-representation "
            "robustness rule: its joint morphology is unusual in the >=2-of-3 shifted-pair "
            "null and both events retain the required positive-node coherence in at least two "
            "representations.  This supports processing robustness but is not yet Level-1 "
            "or physical shock validation."
        )
    else:
        verdict = (
            "The ordered #9 -> #12 pair does not pass the complete frozen "
            "alternative-representation robustness rule.  The earlier BFF pair result is "
            "therefore sensitive to image representation or lacks the required multi-node "
            "coherence in one of the events."
        )

    summary = {
        "events": {
            name: {
                "event_number": candidate["event_number"],
                "utc": candidate["utc"],
                "model": candidate["model"],
                "nodes_rsun": candidate["nodes"].tolist(),
            }
            for name, candidate in CANDIDATES.items()
        },
        "pair_separation_min": PAIR_SEPARATION_MIN,
        "representations": representation_summary,
        "valid_pair_count": int(len(valid_rows)),
        "primary_pair_exceedances": primary_exceedances,
        "primary_pair_p_raw": primary_p,
        "event9_positive_representation_count": event9_positive_reps,
        "event12_positive_representation_count": event12_positive_reps,
        "passed_frozen_robustness_rule": passed,
        "verdict": verdict,
        "no_additional_time_or_pa_shift": True,
        "no_BH": True,
        "data_status": "alternative transformations of Level-0-derived fitpol sector data",
    }

    maps_figure = make_maps_figure(event_results, radii)
    null_figure = make_null_figure(valid_rows, event_results, summary)
    report = write_report(summary)
    summary["files"] = {
        "maps_figure": maps_figure.name,
        "null_figure": null_figure.name,
        "report": report.name,
    }
    json_path = OUT / "pds_20080113_two_event_representation_robustness_results.json"
    json_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    csv_path = OUT / "pds_20080113_two_event_representation_pair_null.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(pair_rows[0]))
        writer.writeheader()
        writer.writerows(pair_rows)
    maps_path = OUT / "pds_20080113_two_event_representation_maps.npz"
    payload = {"radii": radii, "pa_offsets": x2d.PA_OFFSETS}
    for candidate_name in CANDIDATES:
        for representation in REPRESENTATIONS:
            for product in ("total_b", "pb"):
                payload[f"{candidate_name}_{representation}_{product}"] = (
                    event_results[candidate_name][representation]["maps"][product]
                )
    np.savez_compressed(maps_path, **payload)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
