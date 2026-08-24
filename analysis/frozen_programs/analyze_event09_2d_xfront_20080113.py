#!/usr/bin/env python3
"""Frozen 2-D opposite-slope/X-front test for COR1-A event #9.

The test specification is stored in
PDS_20080113_event09_2D_X_test_preregistered.md.  Event #9, its propagation
model, four candidate axial nodes, PA span, line slopes, vertex tolerance and
sideband separation are not fitted after viewing the 2-D event map.
"""

from __future__ import annotations

import csv
import json
import math
from datetime import timedelta
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import gaussian_filter1d

import analyze_event12_shock_cells_20080111_14 as cells
import analyze_pds_event_phase_jitter_20080111_14 as phase
import analyze_pds_nonlinear_transport_20080111_14 as core


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "results" / "pds_20080111_14_nonlinear_level0"
OUT = SOURCE / "event09_2d_xfront"
OUT.mkdir(parents=True, exist_ok=True)

CUBE_FILE = SOURCE / "cor1a_level0_fitpol_sector_15min.npz"
MAPS_FILE = SOURCE / "pds_20080111_14_level0_nonlinear_transport_maps.npz"
PHASE_RESULTS = SOURCE / "event_phase_jitter" / "pds_20080111_14_event_phase_jitter_results.json"

EVENT_NUMBER = 9
EVENT_INDEX = EVENT_NUMBER - 1
EVENT_UTC = core.STRONG_EVENTS[EVENT_INDEX][0]
FIXED_NODES = np.array([1.775, 2.125, 2.500, 2.825], float)
FIXED_MODEL = {
    "family": "deceleration",
    "v_inner_km_s": 100.0,
    "v_outer_km_s": 30.0,
    "change_radius_rsun": 2.20,
}

PA_OFFSETS = np.arange(-12.0, 12.1, 1.0)
ARM_MASK = (np.abs(PA_OFFSETS) >= 2.0) & (np.abs(PA_OFFSETS) <= 8.0)
SLOPES = np.array([0.0250, 0.0375, 0.0500, 0.0625, 0.0750], float)
VERTEX_OFFSETS = np.arange(-0.075, 0.0751, 0.025)
SIDEBAND_RSUN = 0.10
MIN_ARM_SAMPLES = 8


def robust_bff_cube(cube, cadence_min=15.0):
    x = np.asarray(cube, float)
    short_sigma = max(75.0 / cadence_min / 2.355, 0.75)
    long_sigma = max(618.0 / cadence_min / 2.355, 1.5)
    residual = gaussian_filter1d(x, short_sigma, axis=0, mode="nearest")
    residual -= gaussian_filter1d(x, long_sigma, axis=0, mode="nearest")
    center = np.nanmedian(residual, axis=0)
    mad = 1.4826 * np.nanmedian(np.abs(residual - center[None]), axis=0)
    std = np.nanstd(residual, axis=0)
    scale = np.where(mad > 1e-8, mad, np.where(std > 1e-8, std, 1.0))
    return np.clip((residual - center[None]) / scale[None], -6.0, 6.0)


def path_offset_cube(zcube, radii_full, pa, path_full, radii_work):
    radial_indices = np.array([int(np.argmin(np.abs(radii_full - r))) for r in radii_work])
    out = np.empty((zcube.shape[0], len(radii_work), len(PA_OFFSETS)), float)
    pa_step = float(np.median(np.diff(pa)))
    for j, ridx in enumerate(radial_indices):
        targets = path_full[ridx] + PA_OFFSETS
        position = (targets - pa[0]) / pa_step
        lo = np.clip(np.floor(position).astype(int), 0, len(pa) - 1)
        hi = np.clip(lo + 1, 0, len(pa) - 1)
        weight = position - np.floor(position)
        out[:, j, :] = (
            (1.0 - weight[None]) * zcube[:, ridx, lo]
            + weight[None] * zcube[:, ridx, hi]
        )
        outside = (targets < pa[0]) | (targets > pa[-1])
        out[:, j, outside] = np.nan
    return out


def aligned_slice(offset_cube, grid, anchor, tau):
    target = anchor + np.asarray(tau, float)
    dt = float(np.median(np.diff(grid)))
    position = (target - grid[0]) / dt
    valid = (position >= 0) & (position <= len(grid) - 1)
    lo_raw = np.floor(position)
    weight = position - lo_raw
    lo = np.clip(lo_raw.astype(int), 0, len(grid) - 1)
    hi = np.clip(lo + 1, 0, len(grid) - 1)
    ridx = np.arange(len(tau))
    result = (
        (1.0 - weight[:, None]) * offset_cube[lo, ridx, :]
        + weight[:, None] * offset_cube[hi, ridx, :]
    )
    result[~valid] = np.nan
    return result


def sample_map(map2d, radii, radius_targets, pa_indices):
    values = np.empty(len(radius_targets), float)
    for k, (radius, pidx) in enumerate(zip(radius_targets, pa_indices)):
        values[k] = np.interp(radius, radii, map2d[:, pidx], left=np.nan, right=np.nan)
    return values


def diagonal_response(map2d, radii, vertex, signed_slope):
    pa_indices = np.flatnonzero(ARM_MASK)
    q = PA_OFFSETS[pa_indices]
    target = vertex + signed_slope * q
    line = sample_map(map2d, radii, target, pa_indices)
    side_plus = sample_map(map2d, radii, target + SIDEBAND_RSUN, pa_indices)
    side_minus = sample_map(map2d, radii, target - SIDEBAND_RSUN, pa_indices)
    valid = np.isfinite(line) & np.isfinite(side_plus) & np.isfinite(side_minus)
    if np.sum(valid) < MIN_ARM_SAMPLES:
        return np.nan
    contrast = line[valid] - 0.5 * (side_plus[valid] + side_minus[valid])
    return float(np.mean(contrast))


def node_x_response(maps, radii, nominal_node):
    best = None
    for vertex in nominal_node + VERTEX_OFFSETS:
        for slope in SLOPES:
            product_x = {}
            diagonal = {}
            for product, map2d in maps.items():
                plus = diagonal_response(map2d, radii, vertex, slope)
                minus = diagonal_response(map2d, radii, vertex, -slope)
                diagonal[product] = {"positive_slope": plus, "negative_slope": minus}
                product_x[product] = min(plus, minus) if np.isfinite(plus) and np.isfinite(minus) else np.nan
            coherent = min(product_x.values()) if all(np.isfinite(list(product_x.values()))) else np.nan
            candidate = {
                "nominal_radius_rsun": float(nominal_node),
                "vertex_radius_rsun": float(vertex),
                "absolute_slope_rsun_per_degree": float(slope),
                "tb_x_response": float(product_x["total_b"]),
                "pb_x_response": float(product_x["pb"]),
                "coherent_x_response": float(coherent),
                "diagonal_responses": diagonal,
            }
            if best is None or (
                np.isfinite(candidate["coherent_x_response"])
                and (
                    not np.isfinite(best["coherent_x_response"])
                    or candidate["coherent_x_response"] > best["coherent_x_response"]
                )
            ):
                best = candidate
    return best


def score_nodes(maps, radii, nodes, select_top_four=False):
    records = [node_x_response(maps, radii, node) for node in nodes]
    values = np.array([row["coherent_x_response"] for row in records], float)
    finite = values[np.isfinite(values)]
    if select_top_four and len(finite) >= 4:
        used = np.sort(finite)[-4:]
    elif not select_top_four and len(finite) == len(nodes):
        used = finite
    else:
        used = np.array([], float)
    median = float(np.median(used)) if len(used) else np.nan
    return {
        "nodes": records,
        "median_coherent_x_response": median,
        "positive_node_count": int(np.sum(values > 0)),
        "used_node_count": int(len(used)),
    }


def empirical_p(observed, null):
    values = np.asarray(null, float)
    values = values[np.isfinite(values)]
    return float((np.sum(values >= observed) + 1) / (len(values) + 1))


def anchor_utc(anchor):
    return (core.T0 + timedelta(minutes=float(anchor))).strftime("%Y-%m-%dT%H:%M:%SZ")


def figure(event_maps, radii, event_fixed, fixed_null, event_pipeline, joint_p, scan_p):
    fig, axes = plt.subplots(2, 2, figsize=(14, 11), constrained_layout=True)
    products = (("total_b", "tB"), ("pb", "pB"))
    for ax, (product, label) in zip(axes[0], products):
        image = ax.pcolormesh(
            PA_OFFSETS, radii, event_maps[product], cmap="RdBu_r",
            vmin=-2.5, vmax=2.5, shading="auto",
        )
        ax.axvline(0, color="0.2", lw=1.0, ls=":")
        for node in event_fixed["nodes"]:
            q = np.linspace(-8, 8, 81)
            vertex = node["vertex_radius_rsun"]
            slope = node["absolute_slope_rsun_per_degree"]
            ax.plot(q, vertex + slope * q, color="#f6bd60", lw=1.3)
            ax.plot(q, vertex - slope * q, color="#84a59d", lw=1.3)
            ax.scatter([0], [node["nominal_radius_rsun"]], s=24, c="black", zorder=4)
        ax.set(
            title=f"Event #9 travel-aligned {label}",
            xlabel="PA offset from traced ray (deg)", ylabel=r"Radius ($R_\odot$)",
        )
        fig.colorbar(image, ax=ax, label="Temporal BFF robust z")

    same_sign = event_maps["total_b"] * event_maps["pb"] > 0
    coherent_signed = np.zeros_like(event_maps["total_b"])
    coherent_signed[same_sign] = np.sign(event_maps["total_b"][same_sign]) * np.sqrt(
        np.abs(event_maps["total_b"][same_sign] * event_maps["pb"][same_sign])
    )
    image = axes[1, 0].pcolormesh(
        PA_OFFSETS, radii, coherent_signed, cmap="RdBu_r", vmin=-2.5, vmax=2.5,
        shading="auto",
    )
    for node in event_fixed["nodes"]:
        axes[1, 0].axhline(node["nominal_radius_rsun"], color="black", lw=0.8, ls=":")
    axes[1, 0].set(
        title="Same-sign tB/pB geometric-mean map",
        xlabel="PA offset from traced ray (deg)", ylabel=r"Radius ($R_\odot$)",
    )
    fig.colorbar(image, ax=axes[1, 0], label="Coherent signed z")

    finite_null = np.asarray(fixed_null, float)
    finite_null = finite_null[np.isfinite(finite_null)]
    axes[1, 1].hist(finite_null, bins=18, color="0.72", edgecolor="white")
    axes[1, 1].axvline(
        event_fixed["median_coherent_x_response"], color="#9b2c65", lw=2.5,
        label=f"event #9 fixed median={event_fixed['median_coherent_x_response']:.3f}",
    )
    axes[1, 1].set(
        title=(f"Fixed-template null p={event_fixed['p_raw']:.3f}; "
               f"joint p={joint_p:.3f}, 12-anchor={scan_p:.3f}"),
        xlabel="Median coherent X response", ylabel="Unrelated anchors",
    )
    axes[1, 1].legend()
    axes[1, 1].grid(axis="y", alpha=0.2)
    fig.suptitle(
        f"COR1-A event #9 2-D frozen X-front test | pipeline median "
        f"{event_pipeline['median_coherent_x_response']:.3f}", fontsize=15,
    )
    path = OUT / "pds_20080113_1537_event09_2d_xfront_diagnostic.png"
    fig.savefig(path, dpi=230)
    plt.close(fig)
    return path


def write_report(result):
    fixed = result["event9_fixed_template"]
    pipeline = result["event9_pipeline"]
    lines = [
        "# Event #9 frozen 2-D opposite-slope/X-front test",
        "",
        "## Verdict",
        "",
        result["verdict"],
        "",
        "## Frozen-node responses",
        "",
        "| Node (R_sun) | Fitted vertex | |slope| (R_sun/deg) | tB X | pB X | coherent X |",
        "| ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in fixed["nodes"]:
        lines.append(
            f"| {row['nominal_radius_rsun']:.3f} | {row['vertex_radius_rsun']:.3f} | "
            f"{row['absolute_slope_rsun_per_degree']:.4f} | {row['tb_x_response']:.3f} | "
            f"{row['pb_x_response']:.3f} | {row['coherent_x_response']:.3f} |"
        )
    lines += [
        "",
        "## Statistics (raw, no BH)",
        "",
        f"- Fixed-template median coherent X response: {fixed['median_coherent_x_response']:.4f}; "
        f"positive nodes: {fixed['positive_node_count']}/4; p={fixed['p_raw']:.4f} against "
        f"{result['null_anchor_count']} unrelated anchors.",
        f"- Pipeline median of four strongest coherent X responses: "
        f"{pipeline['median_coherent_x_response']:.4f}.",
        f"- Event #9 one-dimensional count-plus-regularity statistic: "
        f"{pipeline['cell_stat']:.4f}.",
        f"- Joint null exceedances: {result['joint_null_exceedance_count']}/"
        f"{result['null_anchor_count']}; raw empirical joint p={result['joint_p_raw']:.4f}.",
        f"- Conservative probability of at least one joint exceedance in a 12-anchor scan: "
        f"{result['joint_scan_probability_12']:.4f}.",
        "",
        "## Meaning",
        "",
        "The statistic tests whether two opposite diagonal line contrasts coexist near each axial "
        "node in both tB and pB.  It does not infer a density jump, magnetosonic Mach number, or "
        "Rankine--Hugoniot consistency.  A visually suggestive polar-map pattern is not called a "
        "shock diamond unless the frozen and joint null tests also pass.",
        "",
        "## Fixed specification",
        "",
        "The event, propagation model, four radii, PA offsets, slope grid, vertex tolerance, "
        "sideband separation, null construction, and interpretation rule were stored in "
        "PDS_20080113_event09_2D_X_test_preregistered.md before the event map was rendered.",
        "",
        "## Calibration limitation",
        "",
        "Inputs are direct bias/exposure-normalized Level-0 tangential polarization fits rather "
        "than complete SECCHI_PREP Level-1 products.  This is a morphology/timing diagnostic.",
    ]
    path = OUT / "PDS_20080113_event09_2D_Xfront_test_no_BH.md"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def main():
    with np.load(CUBE_FILE) as data:
        cube_data = {name: np.asarray(data[name], float) for name in ("total_b", "pb")}
        minutes = np.asarray(data["minutes"], float)
        radii_full = np.asarray(data["radii"], float)
        pa = np.asarray(data["pa"], float)
        quality = np.asarray(data["fit_quality"], float)
    with np.load(MAPS_FILE) as data:
        grid_1d = np.asarray(data["minutes"], float)
        radii = np.asarray(data["radii"], float)
        path_full = np.asarray(data["traced_path_pa"], float)
        zmaps_1d = {
            "total_b": np.asarray(data["z_total_b"], float),
            "pb": np.asarray(data["z_pb"], float),
        }
        raw_1d = {
            "total_b": np.asarray(data["raw_total_b"], float),
            "pb": np.asarray(data["raw_pb"], float),
        }
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
    for product, cube in cube_data.items():
        flat = cube.reshape(len(cube), -1)
        interpolated = np.empty((len(grid), flat.shape[1]), float)
        for j in range(flat.shape[1]):
            valid = good & np.isfinite(flat[:, j])
            interpolated[:, j] = (
                np.interp(grid, minutes[valid], flat[valid, j]) if np.sum(valid) >= 2 else np.nan
            )
        regular[product] = interpolated.reshape((len(grid),) + cube.shape[1:])
        del interpolated

    offset_cubes = {}
    for product, cube in regular.items():
        zcube = robust_bff_cube(cube, core.CADENCE_MIN)
        offset_cubes[product] = path_offset_cube(
            zcube, radii_full, pa, path_full, radii
        )
        del zcube

    fixed_model = phase.model_from_dict(FIXED_MODEL, radii)
    event_anchor = float(event_times[EVENT_INDEX])
    event_maps_fixed = {
        product: aligned_slice(cube, grid, event_anchor, fixed_model.tau)
        for product, cube in offset_cubes.items()
    }
    event_fixed = score_nodes(event_maps_fixed, radii, FIXED_NODES)

    raw_zmaps_1d = {name: cells.robust_standardize(raw) for name, raw in raw_1d.items()}
    families = core.model_families(radii)
    models = []
    for family in ("constant", "acceleration", "deceleration"):
        models.extend(families[family])
    tau_matrix = np.stack([model.tau for model in models])
    phase_results = json.loads(PHASE_RESULTS.read_text())
    reference = phase.model_from_dict(phase_results["reference_model"], radii)
    margin = abs(phase.PHASE_GRID[0]) + abs(phase.LOCAL_GRID[0])
    lo = grid_1d[0] - np.min(reference.tau) + margin
    hi = grid_1d[-1] - margin
    field_anchors = np.arange(
        math.ceil(lo / 30.0) * 30.0,
        math.floor(hi / 30.0) * 30.0 + 0.1,
        30.0,
    )

    event_pipeline_info = cells.analyze_anchor(
        event_anchor, grid_1d, radii, zmaps_1d, raw_zmaps_1d, models, tau_matrix
    )
    event_pipeline_model = phase.model_from_dict(event_pipeline_info["best_model"], radii)
    event_maps_pipeline = {
        product: aligned_slice(cube, grid, event_anchor, event_pipeline_model.tau)
        for product, cube in offset_cubes.items()
    }
    event_nodes_pipeline = np.asarray(
        event_pipeline_info["configs"]["strict"]["matched_primary_radii_rsun"], float
    )
    event_pipeline = score_nodes(
        event_maps_pipeline, radii, event_nodes_pipeline, select_top_four=True
    )
    event_pipeline["cell_stat"] = float(
        event_pipeline_info["configs"]["strict"]["cell_stat"]
    )

    null_rows = []
    for number, anchor in enumerate(field_anchors, 1):
        fixed_maps = {
            product: aligned_slice(cube, grid, float(anchor), fixed_model.tau)
            for product, cube in offset_cubes.items()
        }
        fixed_score = score_nodes(fixed_maps, radii, FIXED_NODES)

        pipeline_info = cells.analyze_anchor(
            float(anchor), grid_1d, radii, zmaps_1d, raw_zmaps_1d, models, tau_matrix
        )
        pipeline_model = phase.model_from_dict(pipeline_info["best_model"], radii)
        pipeline_maps = {
            product: aligned_slice(cube, grid, float(anchor), pipeline_model.tau)
            for product, cube in offset_cubes.items()
        }
        pipeline_nodes = np.asarray(
            pipeline_info["configs"]["strict"]["matched_primary_radii_rsun"], float
        )
        pipeline_score = score_nodes(
            pipeline_maps, radii, pipeline_nodes, select_top_four=True
        )
        null_rows.append({
            "anchor_min": float(anchor),
            "utc": anchor_utc(anchor),
            "fixed_median_x": fixed_score["median_coherent_x_response"],
            "fixed_positive_nodes": fixed_score["positive_node_count"],
            "pipeline_cell_stat": float(pipeline_info["configs"]["strict"]["cell_stat"]),
            "pipeline_matched_peak_count": int(
                pipeline_info["configs"]["strict"]["matched_peak_count"]
            ),
            "pipeline_median_top4_x": pipeline_score["median_coherent_x_response"],
            "best_product": pipeline_info["best_product"],
            "best_model": pipeline_info["best_model"],
        })
        if number % 25 == 0:
            print(f"2-D null: {number}/{len(field_anchors)}", flush=True)

    fixed_null = np.array([row["fixed_median_x"] for row in null_rows], float)
    event_fixed["p_raw"] = empirical_p(event_fixed["median_coherent_x_response"], fixed_null)

    joint_flags = np.array([
        row["pipeline_cell_stat"] >= event_pipeline["cell_stat"]
        and np.isfinite(row["pipeline_median_top4_x"])
        and row["pipeline_median_top4_x"] >= event_pipeline["median_coherent_x_response"]
        for row in null_rows
    ], bool)
    joint_k = int(np.sum(joint_flags))
    joint_p = float((joint_k + 1) / (len(null_rows) + 1))
    scan_probability = float(1.0 - (1.0 - joint_p) ** 12)
    for row, flag in zip(null_rows, joint_flags):
        row["joint_exceedance"] = bool(flag)

    positive_requirement = event_fixed["positive_node_count"] >= 3
    passed = event_fixed["p_raw"] < 0.05 and joint_p < 0.05 and positive_requirement
    if passed:
        verdict = (
            "Event #9 passes the preregistered 2-D morphology criteria: coherent opposite-slope "
            "responses occur at multiple frozen nodes and are unusual in both fixed-template and "
            "joint pipeline nulls.  This supports a shock-cell-like morphology, but Level-1 "
            "calibration and a Mach-number/jump test are still required for a physical shock claim."
        )
    else:
        verdict = (
            "Event #9 does not pass the complete preregistered 2-D detection rule.  Any visible "
            "opposite-slope structure remains a morphology candidate, not a shock-diamond "
            "detection, because the fixed-template rarity, joint pipeline rarity, or multi-node "
            "coherence requirement fails."
        )

    result = {
        "event": "2008-01-13T15:37:30Z",
        "event_number": EVENT_NUMBER,
        "fixed_model": FIXED_MODEL,
        "fixed_nodes_rsun": FIXED_NODES.tolist(),
        "test_parameters": {
            "pa_offsets_deg": PA_OFFSETS.tolist(),
            "arm_pa_absolute_range_deg": [2.0, 8.0],
            "absolute_slopes_rsun_per_degree": SLOPES.tolist(),
            "vertex_offsets_rsun": VERTEX_OFFSETS.tolist(),
            "sideband_rsun": SIDEBAND_RSUN,
            "minimum_arm_samples": MIN_ARM_SAMPLES,
        },
        "event9_fixed_template": event_fixed,
        "event9_pipeline": event_pipeline,
        "null_anchor_count": int(len(null_rows)),
        "joint_null_exceedance_count": joint_k,
        "joint_p_raw": joint_p,
        "joint_scan_probability_12": scan_probability,
        "passed_preregistered_detection_rule": bool(passed),
        "verdict": verdict,
        "no_BH": True,
        "calibration": "bias/exposure-normalized Level-0 tangential fitpol diagnostic",
    }

    figure_path = figure(
        event_maps_fixed, radii, event_fixed, fixed_null, event_pipeline, joint_p,
        scan_probability,
    )
    report_path = write_report(result)
    result["files"] = {"figure": figure_path.name, "report": report_path.name}
    json_path = OUT / "pds_20080113_event09_2d_xfront_results_no_BH.json"
    json_path.write_text(json.dumps(result, indent=2), encoding="utf-8")

    csv_path = OUT / "pds_20080113_event09_2d_xfront_null_table.csv"
    flat_rows = []
    for row in null_rows:
        model = row["best_model"]
        flat_rows.append({
            **{k: v for k, v in row.items() if k != "best_model"},
            "best_family": model["family"],
            "best_v_inner_km_s": model["v_inner_km_s"],
            "best_v_outer_km_s": model["v_outer_km_s"],
            "best_change_radius_rsun": model["change_radius_rsun"],
        })
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(flat_rows[0]))
        writer.writeheader()
        writer.writerows(flat_rows)

    maps_path = OUT / "pds_20080113_1537_event09_2d_xfront_maps.npz"
    np.savez_compressed(
        maps_path, radii=radii, pa_offsets=PA_OFFSETS,
        event_total_b=event_maps_fixed["total_b"], event_pb=event_maps_fixed["pb"],
        fixed_null_median_x=fixed_null,
    )
    print(json.dumps({
        "event_fixed_median_x": event_fixed["median_coherent_x_response"],
        "event_fixed_positive_nodes": event_fixed["positive_node_count"],
        "event_fixed_p_raw": event_fixed["p_raw"],
        "event_pipeline_cell_stat": event_pipeline["cell_stat"],
        "event_pipeline_median_top4_x": event_pipeline["median_coherent_x_response"],
        "joint_null_exceedances": joint_k,
        "joint_p_raw": joint_p,
        "joint_scan_probability_12": scan_probability,
        "passed": passed,
        "node_responses": event_fixed["nodes"],
    }, indent=2))


if __name__ == "__main__":
    main()
