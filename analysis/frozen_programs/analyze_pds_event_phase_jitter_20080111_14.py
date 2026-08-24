#!/usr/bin/env python3
"""Event-by-event phase-jitter diagnostic for 11--14 January 2008.

The 12 COR2 anchors, PA, Level-0 COR1 maps, BFF filtering and model grid are
inherited from the frozen nonlinear-transport test.  Raw empirical p values
are reported without BH adjustment.

Two complementary questions are asked.

1. Does each event contain an identifiable COR1 ridge after the complete
   product/kinematic-model scan is repeated at unrelated anchor times?
2. Relative to one common nonlinear reference path, does each event require a
   different global phase offset or a height-dependent phase kink?

The phase-kink statistic is descriptive unless the event ridge itself is
identified.  Local maxima in noise can otherwise manufacture an apparent
phase curve.
"""

from __future__ import annotations

import csv
import json
import math
from dataclasses import asdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import analyze_pds_nonlinear_transport_20080111_14 as core


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "results" / "pds_20080111_14_nonlinear_level0"
OUT = SOURCE / "event_phase_jitter"
OUT.mkdir(parents=True, exist_ok=True)
MAPS_FILE = SOURCE / "pds_20080111_14_level0_nonlinear_transport_maps.npz"
BASE_RESULTS_FILE = SOURCE / "pds_20080111_14_level0_nonlinear_transport_results.json"

PHASE_GRID = np.arange(-90.0, 90.1, core.CADENCE_MIN)
LOCAL_GRID = np.arange(-60.0, 60.1, core.CADENCE_MIN)
KINK_RADII = np.arange(1.8, 2.901, 0.05)


def empirical_p(observed, null):
    values = np.asarray(null, float)
    values = values[np.isfinite(values)]
    return float((np.sum(values >= observed) + 1) / (len(values) + 1))


def model_from_dict(row, radii):
    family = row["family"]
    vin = float(row["v_inner_km_s"])
    vout = float(row["v_outer_km_s"])
    rc = row["change_radius_rsun"]
    if rc is None:
        tau = -(3.0 - radii) * core.RSUN_KM / vin / 60.0
    else:
        tau = core.travel_time_piecewise(radii, vin, vout, float(rc))
    return core.Model(family, vin, vout, None if rc is None else float(rc), tau)


def interpolated_ridge(zmap, grid, anchor, tau):
    target = anchor + tau
    values = np.array([
        np.interp(target[j], grid, zmap[:, j], left=np.nan, right=np.nan)
        for j in range(len(tau))
    ])
    return values


def score_one_model(zmap, grid, anchor, tau):
    values = interpolated_ridge(zmap, grid, anchor, tau)
    return float(np.nanmean(values)) if np.mean(np.isfinite(values)) >= 0.80 else -np.inf


def vector_model_scores(zmap, grid, anchor, tau_matrix):
    dt = float(np.median(np.diff(grid)))
    pos = (anchor + tau_matrix - grid[0]) / dt
    valid = (pos >= 0.0) & (pos <= len(grid) - 1)
    lo_raw = np.floor(pos)
    weight = pos - lo_raw
    lo = np.clip(lo_raw.astype(int), 0, len(grid) - 1)
    hi = np.clip(lo + 1, 0, len(grid) - 1)
    ridx = np.arange(zmap.shape[1])[None, :]
    values = (1.0 - weight) * zmap[lo, ridx] + weight * zmap[hi, ridx]
    values = np.where(valid, values, np.nan)
    with np.errstate(invalid="ignore"):
        scores = np.nanmean(values, axis=1)
    scores[np.mean(valid, axis=1) < 0.80] = -np.inf
    return scores


def best_scanned_model(zmaps, grid, anchor, models, tau_matrix):
    winner = None
    for product, zmap in zmaps.items():
        scores = vector_model_scores(zmap, grid, anchor, tau_matrix)
        index = int(np.nanargmax(scores))
        candidate = (float(scores[index]), product, index, models[index])
        if winner is None or candidate[0] > winner[0]:
            winner = candidate
    return winner


def best_global_phase(zmap, grid, anchor, tau):
    scores = np.array([
        score_one_model(zmap, grid, anchor + delta, tau) for delta in PHASE_GRID
    ])
    index = int(np.nanargmax(scores))
    return float(PHASE_GRID[index]), float(scores[index]), scores


def local_peak_phase_curve(zmap, grid, anchor, tau):
    offsets = np.empty(len(tau), float)
    peaks = np.empty(len(tau), float)
    for j in range(len(tau)):
        target = anchor + tau[j] + LOCAL_GRID
        values = np.interp(target, grid, zmap[:, j], left=np.nan, right=np.nan)
        if not np.any(np.isfinite(values)):
            offsets[j] = np.nan
            peaks[j] = np.nan
            continue
        k = int(np.nanargmax(values))
        offset = float(LOCAL_GRID[k])
        peak = float(values[k])
        if 0 < k < len(values) - 1 and np.all(np.isfinite(values[k-1:k+2])):
            ym, y0, yp = values[k-1:k+2]
            curvature = ym - 2.0 * y0 + yp
            if curvature < -1e-8:
                fraction = float(np.clip(0.5 * (ym - yp) / curvature, -0.5, 0.5))
                offset += fraction * core.CADENCE_MIN
                peak = float(y0 - 0.25 * (ym - yp) * fraction)
        offsets[j] = offset
        peaks[j] = peak
    return offsets, peaks


def weighted_mean(values, weights):
    valid = np.isfinite(values) & np.isfinite(weights) & (weights > 0)
    if not np.any(valid):
        return np.nan
    return float(np.sum(values[valid] * weights[valid]) / np.sum(weights[valid]))


def weighted_kink(radii, offsets, peaks):
    valid = np.isfinite(offsets) & np.isfinite(peaks)
    weights = np.square(np.clip(peaks, 0.0, 3.0))
    if np.sum(valid & (weights > 0)) < 12:
        weights = np.where(valid, 1.0, 0.0)
    center = weighted_mean(offsets, weights)
    if not np.isfinite(center):
        return {"radius": np.nan, "jump": np.nan, "gain": np.nan,
                "inner_phase": np.nan, "outer_phase": np.nan}
    sse0 = float(np.sum(weights[valid] * (offsets[valid] - center) ** 2))
    best = None
    for rc in KINK_RADII:
        inner = valid & (radii < rc) & (weights > 0)
        outer = valid & (radii >= rc) & (weights > 0)
        if np.sum(inner) < 6 or np.sum(outer) < 6:
            continue
        phase_inner = weighted_mean(offsets[inner], weights[inner])
        phase_outer = weighted_mean(offsets[outer], weights[outer])
        sse = float(np.sum(weights[inner] * (offsets[inner] - phase_inner) ** 2)
                    + np.sum(weights[outer] * (offsets[outer] - phase_outer) ** 2))
        candidate = (sse, float(rc), float(phase_outer - phase_inner),
                     float(phase_inner), float(phase_outer))
        if best is None or candidate[0] < best[0]:
            best = candidate
    if best is None:
        return {"radius": np.nan, "jump": np.nan, "gain": np.nan,
                "inner_phase": np.nan, "outer_phase": np.nan}
    gain = max(0.0, (sse0 - best[0]) / sse0) if sse0 > 1e-8 else 0.0
    return {"radius": best[1], "jump": best[2], "gain": float(gain),
            "inner_phase": best[3], "outer_phase": best[4]}


def robust_mad(values):
    values = np.asarray(values, float)
    values = values[np.isfinite(values)]
    if len(values) == 0:
        return np.nan
    center = np.median(values)
    return float(1.4826 * np.median(np.abs(values - center)))


def model_record(model):
    return {
        "family": model.family,
        "v_inner_km_s": float(model.v_inner),
        "v_outer_km_s": float(model.v_outer),
        "change_radius_rsun": None if model.change_radius is None else float(model.change_radius),
    }


def make_phase_curve_figure(rows, radii, curves_t, peaks_t, curves_p):
    fig, axes = plt.subplots(3, 4, figsize=(16, 11.5), sharex=True, sharey=True,
                             constrained_layout=True)
    for i, ax in enumerate(axes.flat):
        row = rows[i]
        strong = peaks_t[i] >= 0.75
        ax.axhline(0.0, color="0.55", lw=0.8)
        ax.axvline(2.10, color="#e09f3e", ls="--", lw=1.2, label="common 2.10")
        ax.plot(radii, curves_p[i], color="#4c78a8", lw=1.0, ls="--", alpha=0.65,
                label="pB residual")
        ax.plot(radii, curves_t[i], color="#a23b72", lw=1.5, label="tB residual")
        ax.scatter(radii[~strong], curves_t[i][~strong], s=9, color="0.72", alpha=0.55)
        ax.scatter(radii[strong], curves_t[i][strong], s=13, color="#a23b72")
        if np.isfinite(row["kink_radius_rsun"]):
            ax.axvline(row["kink_radius_rsun"], color="#2a9d8f", ls=":", lw=1.2)
        ax.set_title(
            f"#{i+1:02d} {row['utc'][5:16].replace('T', ' ')}  "
            f"p(ridge)={row['best_ridge_p_raw']:.3f}\n"
            f"kink={row['kink_radius_rsun']:.2f} R$_\\odot$, p={row['kink_p_raw']:.3f}",
            fontsize=9,
        )
        ax.set_ylim(-68, 68)
        ax.grid(alpha=0.15)
    for ax in axes[-1]:
        ax.set_xlabel(r"Radius ($R_\odot$)")
    for ax in axes[:, 0]:
        ax.set_ylabel("Residual phase (min)")
    axes[0, 0].legend(loc="upper right", fontsize=7, frameon=True)
    path = OUT / "pds_20080111_14_event_phase_curves.png"
    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return path


def make_summary_figure(rows, radii, curves_t, global_tests):
    fig, axes = plt.subplots(3, 1, figsize=(13.5, 10.5), constrained_layout=True)
    image = np.stack(curves_t)
    pcm = axes[0].pcolormesh(radii, np.arange(1, 13), image, shading="auto",
                             cmap="RdBu_r", vmin=-60, vmax=60)
    axes[0].axvline(2.10, color="k", ls="--", lw=1)
    axes[0].set(ylabel="Event number", title="Height-dependent tB phase residual after event-global offset")
    axes[0].set_yticks(np.arange(1, 13))
    fig.colorbar(pcm, ax=axes[0], label="Residual phase (min)")

    x = np.arange(1, 13)
    dt = np.array([row["global_phase_tb_min"] for row in rows])
    dp = np.array([row["global_phase_pb_min"] for row in rows])
    axes[1].plot(x, dt, "o-", color="#a23b72", label="tB")
    axes[1].plot(x, dp, "s--", color="#4c78a8", label="pB")
    axes[1].axhline(0.0, color="0.5", lw=0.8)
    axes[1].set(ylabel="Event-global phase offset (min)", title=(
        f"Per-event phase shifts; global jitter-gain p(tB)={global_tests['tb']['gain_p_raw']:.3f}, "
        f"p(pB)={global_tests['pb']['gain_p_raw']:.3f}"))
    axes[1].set_xticks(x)
    axes[1].legend()
    axes[1].grid(alpha=0.2)

    p_best = np.array([row["best_ridge_p_raw"] for row in rows])
    p_common = np.array([row["common_ridge_p_raw"] for row in rows])
    p_kink = np.array([row["kink_p_raw"] for row in rows])
    axes[2].plot(x, -np.log10(p_best), "o-", label="best scanned ridge")
    axes[2].plot(x, -np.log10(p_common), "s--", label="common path")
    axes[2].plot(x, -np.log10(p_kink), "^:", label="phase kink")
    axes[2].axhline(-math.log10(0.05), color="k", ls="--", lw=1, label="raw p=0.05")
    axes[2].set(xlabel="Event number", ylabel=r"$-\log_{10}(p_{raw})$",
                title="Individual raw empirical p values (no BH)")
    axes[2].set_xticks(x)
    axes[2].grid(alpha=0.2)
    axes[2].legend(ncol=4, fontsize=9)
    path = OUT / "pds_20080111_14_event_phase_summary.png"
    fig.savefig(path, dpi=220)
    plt.close(fig)
    return path


def write_report(rows, reference, global_tests, field_null_count):
    ridge_count = sum(row["best_ridge_p_raw"] < 0.05 for row in rows)
    common_count = sum(row["common_ridge_p_raw"] < 0.05 for row in rows)
    kink_count = sum(
        row["best_ridge_p_raw"] < 0.05 and row["kink_p_raw"] < 0.05
        for row in rows
    )
    tb_global_boundary = sum(abs(row["global_phase_tb_min"]) >= 90 for row in rows)
    pb_global_boundary = sum(abs(row["global_phase_pb_min"]) >= 90 for row in rows)
    median_local_boundary = float(np.median([
        row["phase_curve_tb_boundary_fraction"] for row in rows
    ]))
    lines = [
        "# Event-dependent phase-jitter test: 11--14 January 2008",
        "",
        "## Result",
        "",
        f"The full product/model scan identifies {ridge_count}/12 individual COR1 ridges at raw p<0.05; "
        f"the frozen common nonlinear path gives a nominal raw p<0.05 for {common_count}/12.  Only {kink_count}/12 events "
        "simultaneously have an individually identified ridge and a raw-significant height-dependent "
        "phase kink.  These are unadjusted diagnostic counts, not family-wise discoveries.",
        "",
        f"Allowing every event its own global tB phase shift improves the common-path mean score with "
        f"raw common-shift p={global_tests['tb']['gain_p_raw']:.4f}; for pB the corresponding p is "
        f"{global_tests['pb']['gain_p_raw']:.4f}.  The observed robust event-to-event phase scatter is "
        f"{global_tests['tb']['observed_phase_scatter_min']:.1f} min in tB and "
        f"{global_tests['pb']['observed_phase_scatter_min']:.1f} min in pB.",
        "",
        f"The phase optimizer lands on the +/-90 min search boundary for {tb_global_boundary}/12 tB "
        f"events and {pb_global_boundary}/12 pB events; the median fraction of local tB phase points "
        f"at the +/-60 min boundary is {100*median_local_boundary:.1f}%.  This is a direct warning that "
        "many individual phases are not constrained by a ridge.",
        "",
        "The only selection-controlled raw candidate is event #12 (2008-01-13 21:37:30 UT, "
        f"p={rows[11]['best_ridge_p_raw']:.4f}), which prefers a tB deceleration alignment "
        f"{rows[11]['best_model']['v_inner_km_s']:.0f}->{rows[11]['best_model']['v_outer_km_s']:.0f} "
        f"km/s near {rows[11]['best_model']['change_radius_rsun']:.2f} R_sun.  Its phase-kink "
        f"p={rows[11]['kink_p_raw']:.3f}, however, is not significant and both tB and pB global "
        "phases hit the -90 min search boundary.  Event #11, exactly 120 min earlier, is only "
        f"near-significant (p={rows[10]['best_ridge_p_raw']:.4f}) and prefers the opposite family: "
        f"acceleration near {rows[10]['best_model']['change_radius_rsun']:.2f} R_sun.  The pair therefore "
        "does not define a shared moving-shock trajectory; it is also not independent of the previously "
        "discussed 80--130 min variability band.",
        "",
        "## Frozen reference",
        "",
        f"- Common diagnostic path: {reference.family}, {reference.v_inner:.0f}->{reference.v_outer:.0f} "
        f"km/s at {reference.change_radius:.2f} R_sun, selected previously only for visualization.",
        "- 12 COR2 anchors and PA=174.5 deg were fixed independently of this phase analysis.",
        "- Phase search uses 15-min cadence, event-global offsets of +/-90 min and local residual offsets "
        "of +/-60 min.",
        f"- Individual ridge and kink nulls use {field_null_count} unrelated anchor times and repeat all "
        "allowed choices.",
        "- Raw empirical p values are reported; no BH adjustment.",
        "",
        "## Event table",
        "",
        "| # | COR2 anchor UTC | z | Best ridge | p raw | Common-path p | tB/pB global phase (min) | "
        "Kink radius | Jump (min) | Kink p raw | Interpretation |",
        "| ---: | --- | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |",
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
        lines.append(
            f"| {row['event_number']} | {row['utc'].replace('T', ' ')[:16]} | {row['cor2_z']:.2f} | "
            f"{model_text} | {row['best_ridge_p_raw']:.3f} | {row['common_ridge_p_raw']:.3f} | "
            f"{row['global_phase_tb_min']:+.0f}/{row['global_phase_pb_min']:+.0f} | "
            f"{row['kink_radius_rsun']:.2f} | {row['kink_jump_min']:+.1f} | "
            f"{row['kink_p_raw']:.3f} | {row['interpretation']} |"
        )
    lines += [
        "",
        "## Reading the phase curves",
        "",
        "A stationary shock can preserve the incoming cadence while producing the same kink height and "
        "phase/amplitude change in repeated events.  A moving or reforming shock can instead give different "
        "global offsets, kink heights and phase jumps.  However, a phase curve is physically interpretable "
        "only when the corresponding ridge is identified against unrelated times.  Otherwise the local-peak "
        "operator is following maxima in filtered noise.",
        "",
        "The event-specific best models are exploratory.  Their speeds are pattern-alignment speeds and "
        "must not be read as plasma bulk velocities.  A shock claim still requires a repeatable brightness "
        "compression, upstream fast-mode Mach number above unity and Rankine--Hugoniot consistency.",
        "",
        "## Calibration limitation",
        "",
        "The maps use direct bias/exposure-normalized COR1-A Level-0 polarization triplets and not the "
        "complete SECCHI_PREP/background chain.  Absolute pB compression ratios are therefore not inferred.",
    ]
    path = OUT / "PDS_20080111_14_event_phase_jitter_no_BH.md"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def main():
    with np.load(MAPS_FILE) as data:
        grid = np.asarray(data["minutes"], float)
        radii = np.asarray(data["radii"], float)
        zmaps = {
            "total_b": np.asarray(data["z_total_b"], float),
            "pb": np.asarray(data["z_pb"], float),
        }
        event_times = np.asarray(data["event_times"], float)
    base = json.loads(BASE_RESULTS_FILE.read_text(encoding="utf-8"))
    reference = model_from_dict(base["diagnostic_selected_model"], radii)

    core.RADII_WORK = radii
    family_dict = core.model_families(radii)
    models = []
    for family in ("constant", "acceleration", "deceleration"):
        models.extend(family_dict[family])
    tau_matrix = np.stack([model.tau for model in models])

    # Use one field of unrelated times for all event-wise empirical tests.
    phase_margin = abs(PHASE_GRID[0]) + abs(LOCAL_GRID[0])
    lo = grid[0] - np.min(reference.tau) + phase_margin
    hi = grid[-1] - phase_margin
    field_anchors = np.arange(math.ceil(lo / 30.0) * 30.0,
                              math.floor(hi / 30.0) * 30.0 + 0.1, 30.0)
    field_best_scores = []
    field_reference_scores = []
    field_gain_scores = []
    field_kink_gains = []
    for number, anchor in enumerate(field_anchors, 1):
        best = best_scanned_model(zmaps, grid, float(anchor), models, tau_matrix)
        ref = score_one_model(zmaps["total_b"], grid, float(anchor), reference.tau)
        field_best_scores.append(best[0])
        field_reference_scores.append(ref)
        field_gain_scores.append(best[0] - ref)
        delta, _, _ = best_global_phase(zmaps["total_b"], grid, float(anchor), reference.tau)
        offsets, peaks = local_peak_phase_curve(
            zmaps["total_b"], grid, float(anchor) + delta, reference.tau
        )
        field_kink_gains.append(weighted_kink(radii, offsets, peaks)["gain"])
        if number % 50 == 0:
            print(f"Field null: {number}/{len(field_anchors)}", flush=True)
    field_best_scores = np.asarray(field_best_scores)
    field_reference_scores = np.asarray(field_reference_scores)
    field_gain_scores = np.asarray(field_gain_scores)
    field_kink_gains = np.asarray(field_kink_gains)

    rows = []
    curves_t, curves_p, peaks_t_all, peaks_p_all = [], [], [], []
    for i, (utc, cor2_z) in enumerate(core.STRONG_EVENTS):
        anchor = float(event_times[i])
        best_score, product, model_index, best_model = best_scanned_model(
            zmaps, grid, anchor, models, tau_matrix
        )
        common_score = score_one_model(zmaps["total_b"], grid, anchor, reference.tau)
        phase_t, phase_score_t, _ = best_global_phase(
            zmaps["total_b"], grid, anchor, reference.tau
        )
        phase_p, phase_score_p, _ = best_global_phase(
            zmaps["pb"], grid, anchor, reference.tau
        )
        curve_t, peaks_t = local_peak_phase_curve(
            zmaps["total_b"], grid, anchor + phase_t, reference.tau
        )
        curve_p, peaks_p = local_peak_phase_curve(
            zmaps["pb"], grid, anchor + phase_p, reference.tau
        )
        kink = weighted_kink(radii, curve_t, peaks_t)
        ridge_p = empirical_p(best_score, field_best_scores)
        common_p = empirical_p(common_score, field_reference_scores)
        gain_p = empirical_p(best_score - common_score, field_gain_scores)
        kink_p = empirical_p(kink["gain"], field_kink_gains)
        strong = peaks_t >= 0.75
        phase_cross_mad = robust_mad(curve_t - curve_p)
        if ridge_p < 0.05 and kink_p < 0.05:
            interpretation = "ridge + kink candidate"
        elif ridge_p < 0.05:
            interpretation = "ridge; kink not identified"
        else:
            interpretation = "phase not identifiable"
        row = {
            "event_number": i + 1,
            "utc": utc,
            "cor2_anchor_min": anchor,
            "cor2_z": float(cor2_z),
            "best_product": product,
            "best_model": model_record(best_model),
            "best_ridge_score": float(best_score),
            "best_ridge_p_raw": ridge_p,
            "common_ridge_score": float(common_score),
            "common_ridge_p_raw": common_p,
            "individual_model_gain": float(best_score - common_score),
            "individual_model_gain_p_raw": gain_p,
            "global_phase_tb_min": phase_t,
            "global_phase_pb_min": phase_p,
            "global_phase_product_difference_min": abs(phase_t - phase_p),
            "global_phase_tb_score": phase_score_t,
            "global_phase_pb_score": phase_score_p,
            "phase_curve_tb_mad_min": robust_mad(curve_t),
            "phase_curve_tb_pb_difference_mad_min": phase_cross_mad,
            "phase_curve_tb_boundary_fraction": float(np.mean(np.abs(curve_t) >= 59.9)),
            "strong_phase_height_fraction": float(np.mean(strong)),
            "kink_radius_rsun": float(kink["radius"]),
            "kink_jump_min": float(kink["jump"]),
            "kink_gain": float(kink["gain"]),
            "kink_p_raw": kink_p,
            "interpretation": interpretation,
        }
        rows.append(row)
        curves_t.append(curve_t)
        curves_p.append(curve_p)
        peaks_t_all.append(peaks_t)
        peaks_p_all.append(peaks_p)

    curves_t = np.asarray(curves_t)
    curves_p = np.asarray(curves_p)
    peaks_t_all = np.asarray(peaks_t_all)
    peaks_p_all = np.asarray(peaks_p_all)

    # Family-level test: does independent event timing improve the common path
    # more than it does after shifting the complete 12-event train together?
    shift_lo = lo - np.min(event_times)
    shift_hi = hi - np.max(event_times)
    common_shifts = np.arange(math.ceil(shift_lo / core.CADENCE_MIN) * core.CADENCE_MIN,
                              math.floor(shift_hi / core.CADENCE_MIN) * core.CADENCE_MIN + 0.1,
                              core.CADENCE_MIN)
    null_shifts = common_shifts[np.abs(common_shifts) > 1e-6]
    global_tests = {}
    for product, zmap in zmaps.items():
        observed_fixed = []
        observed_best = []
        observed_phases = []
        for anchor in event_times:
            observed_fixed.append(score_one_model(zmap, grid, float(anchor), reference.tau))
            delta, score, _ = best_global_phase(zmap, grid, float(anchor), reference.tau)
            observed_best.append(score)
            observed_phases.append(delta)
        observed_gain = float(np.mean(observed_best) - np.mean(observed_fixed))
        null_gain = []
        null_scatter = []
        for shift in null_shifts:
            fixed, best_values, phases = [], [], []
            for anchor in event_times + shift:
                fixed.append(score_one_model(zmap, grid, float(anchor), reference.tau))
                delta, score, _ = best_global_phase(zmap, grid, float(anchor), reference.tau)
                best_values.append(score)
                phases.append(delta)
            null_gain.append(float(np.mean(best_values) - np.mean(fixed)))
            null_scatter.append(robust_mad(phases))
        key = "tb" if product == "total_b" else "pb"
        global_tests[key] = {
            "observed_fixed_score": float(np.mean(observed_fixed)),
            "observed_event_shifted_score": float(np.mean(observed_best)),
            "observed_gain": observed_gain,
            "gain_p_raw": empirical_p(observed_gain, null_gain),
            "observed_phase_offsets_min": [float(x) for x in observed_phases],
            "observed_phase_scatter_min": robust_mad(observed_phases),
            "phase_scatter_p_raw": empirical_p(robust_mad(observed_phases), null_scatter),
            "null_shift_count": int(len(null_shifts)),
        }

    curve_path = make_phase_curve_figure(rows, radii, curves_t, peaks_t_all, curves_p)
    summary_path = make_summary_figure(rows, radii, curves_t, global_tests)
    report_path = write_report(rows, reference, global_tests, len(field_anchors))

    csv_path = OUT / "pds_20080111_14_event_phase_jitter_table.csv"
    flat_rows = []
    for row in rows:
        flat = {key: value for key, value in row.items() if key != "best_model"}
        flat.update({f"best_{key}": value for key, value in row["best_model"].items()})
        flat_rows.append(flat)
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(flat_rows[0]))
        writer.writeheader()
        writer.writerows(flat_rows)

    npz_path = OUT / "pds_20080111_14_event_phase_jitter_curves.npz"
    np.savez_compressed(
        npz_path,
        radii=radii,
        event_times=event_times,
        event_utc=np.array([item[0] for item in core.STRONG_EVENTS]),
        phase_residual_tb=curves_t.astype(np.float32),
        phase_residual_pb=curves_p.astype(np.float32),
        phase_peak_z_tb=peaks_t_all.astype(np.float32),
        phase_peak_z_pb=peaks_p_all.astype(np.float32),
    )

    output = {
        "event": "2008-01-11--14",
        "reference_model": model_record(reference),
        "event_count": len(rows),
        "field_null_count": int(len(field_anchors)),
        "no_BH": True,
        "global_event_phase_tests": global_tests,
        "events": rows,
        "files": {
            "report": report_path.name,
            "table": csv_path.name,
            "phase_curves": curve_path.name,
            "summary_figure": summary_path.name,
            "curve_data": npz_path.name,
        },
        "calibration_caveat": "Direct Level-0 fitpol diagnostic; no absolute pB compression ratio",
    }
    json_path = OUT / "pds_20080111_14_event_phase_jitter_results.json"
    json_path.write_text(json.dumps(output, indent=2), encoding="utf-8")
    print(json.dumps({
        "global_tests": global_tests,
        "ridge_p_raw": [round(row["best_ridge_p_raw"], 4) for row in rows],
        "common_p_raw": [round(row["common_ridge_p_raw"], 4) for row in rows],
        "kink_p_raw": [round(row["kink_p_raw"], 4) for row in rows],
        "phase_tb": [row["global_phase_tb_min"] for row in rows],
        "phase_pb": [row["global_phase_pb_min"] for row in rows],
    }, indent=2))


if __name__ == "__main__":
    main()
