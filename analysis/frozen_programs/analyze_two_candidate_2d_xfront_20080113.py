#!/usr/bin/env python3
"""Frozen two-candidate 2-D X-front comparison for events #9 and #12."""

from __future__ import annotations

import csv
import json
import math
from datetime import timedelta
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import analyze_event09_2d_xfront_20080113 as x2d
import analyze_event12_shock_cells_20080111_14 as cells
import analyze_pds_event_phase_jitter_20080111_14 as phase
import analyze_pds_nonlinear_transport_20080111_14 as core


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "results" / "pds_20080111_14_nonlinear_level0"
OUT = SOURCE / "two_candidate_2d_xfront"
OUT.mkdir(parents=True, exist_ok=True)

EVENT9_RESULTS = SOURCE / "event09_2d_xfront" / "pds_20080113_event09_2d_xfront_results_no_BH.json"
EVENT9_NULL_TABLE = SOURCE / "event09_2d_xfront" / "pds_20080113_event09_2d_xfront_null_table.csv"
PAIR_SEPARATION_MIN = 360.0

CANDIDATE12 = {
    "event_number": 12,
    "event_index": 11,
    "utc": core.STRONG_EVENTS[11][0],
    "model": {
        "family": "deceleration",
        "v_inner_km_s": 100.0,
        "v_outer_km_s": 25.0,
        "change_radius_rsun": 2.90,
    },
    "nodes_rsun": np.array([1.950, 2.350, 2.700], float),
    "required_positive_nodes": 3,
}


def prepare_offset_cubes():
    with np.load(x2d.CUBE_FILE) as data:
        cubes = {name: np.asarray(data[name], float) for name in ("total_b", "pb")}
        minutes = np.asarray(data["minutes"], float)
        radii_full = np.asarray(data["radii"], float)
        pa = np.asarray(data["pa"], float)
        quality = np.asarray(data["fit_quality"], float)
    with np.load(x2d.MAPS_FILE) as data:
        maps = {
            "grid": np.asarray(data["minutes"], float),
            "radii": np.asarray(data["radii"], float),
            "path_full": np.asarray(data["traced_path_pa"], float),
            "zmaps": {
                "total_b": np.asarray(data["z_total_b"], float),
                "pb": np.asarray(data["z_pb"], float),
            },
            "raw": {
                "total_b": np.asarray(data["raw_total_b"], float),
                "pb": np.asarray(data["raw_pb"], float),
            },
            "event_times": np.asarray(data["event_times"], float),
        }

    center = np.nanmedian(quality)
    mad = 1.4826 * np.nanmedian(np.abs(quality - center))
    good = np.isfinite(quality) & (quality <= center + 6 * max(mad, 1e-8))
    grid = np.arange(
        math.ceil(np.min(minutes[good]) / core.CADENCE_MIN) * core.CADENCE_MIN,
        math.floor(np.max(minutes[good]) / core.CADENCE_MIN) * core.CADENCE_MIN + 0.1,
        core.CADENCE_MIN,
    )
    offsets = {}
    for product, cube in cubes.items():
        flat = cube.reshape(len(cube), -1)
        regular = np.empty((len(grid), flat.shape[1]), float)
        for j in range(flat.shape[1]):
            valid = good & np.isfinite(flat[:, j])
            regular[:, j] = (
                np.interp(grid, minutes[valid], flat[valid, j]) if np.sum(valid) >= 2 else np.nan
            )
        regular = regular.reshape((len(grid),) + cube.shape[1:])
        zcube = x2d.robust_bff_cube(regular, core.CADENCE_MIN)
        offsets[product] = x2d.path_offset_cube(
            zcube, radii_full, pa, maps["path_full"], maps["radii"]
        )
        del regular, zcube
    return grid, offsets, maps


def score_top_n(maps2d, radii, nodes, n_required):
    records = [x2d.node_x_response(maps2d, radii, node) for node in nodes]
    values = np.array([row["coherent_x_response"] for row in records], float)
    finite = values[np.isfinite(values)]
    used = np.sort(finite)[-n_required:] if len(finite) >= n_required else np.array([], float)
    return {
        "nodes": records,
        "median_coherent_x_response": float(np.median(used)) if len(used) else np.nan,
        "positive_node_count": int(np.sum(values > 0)),
        "used_node_count": int(len(used)),
    }


def empirical_p(observed, null):
    values = np.asarray(null, float)
    values = values[np.isfinite(values)]
    return float((np.sum(values >= observed) + 1) / (len(values) + 1))


def anchor_utc(anchor):
    return (core.T0 + timedelta(minutes=float(anchor))).strftime("%Y-%m-%dT%H:%M:%SZ")


def make_figure(event_maps, radii, fixed, fixed_null, comparison):
    fig, axes = plt.subplots(2, 2, figsize=(14, 11), constrained_layout=True)
    for ax, product, label in zip(axes[0], ("total_b", "pb"), ("tB", "pB")):
        image = ax.pcolormesh(
            x2d.PA_OFFSETS, radii, event_maps[product], cmap="RdBu_r",
            vmin=-2.5, vmax=2.5, shading="auto",
        )
        ax.axvline(0, color="0.2", ls=":", lw=1)
        for node in fixed["nodes"]:
            q = np.linspace(-8, 8, 81)
            vertex = node["vertex_radius_rsun"]
            slope = node["absolute_slope_rsun_per_degree"]
            ax.plot(q, vertex + slope * q, color="#f6bd60", lw=1.3)
            ax.plot(q, vertex - slope * q, color="#84a59d", lw=1.3)
            ax.scatter([0], [node["nominal_radius_rsun"]], c="black", s=25, zorder=4)
        ax.set(
            title=f"Event #12 frozen {label}",
            xlabel="PA offset from traced ray (deg)", ylabel=r"Radius ($R_\odot$)",
        )
        fig.colorbar(image, ax=ax, label="Temporal BFF robust z")

    same = event_maps["total_b"] * event_maps["pb"] > 0
    coherent = np.zeros_like(event_maps["total_b"])
    coherent[same] = np.sign(event_maps["total_b"][same]) * np.sqrt(
        np.abs(event_maps["total_b"][same] * event_maps["pb"][same])
    )
    image = axes[1, 0].pcolormesh(
        x2d.PA_OFFSETS, radii, coherent, cmap="RdBu_r", vmin=-2.5, vmax=2.5,
        shading="auto",
    )
    for node in CANDIDATE12["nodes_rsun"]:
        axes[1, 0].axhline(node, color="black", lw=0.8, ls=":")
    axes[1, 0].set(
        title="Event #12 same-sign tB/pB map",
        xlabel="PA offset from traced ray (deg)", ylabel=r"Radius ($R_\odot$)",
    )
    fig.colorbar(image, ax=axes[1, 0], label="Coherent signed z")

    values = np.asarray(fixed_null, float)
    values = values[np.isfinite(values)]
    axes[1, 1].hist(values, bins=18, color="0.72", edgecolor="white")
    axes[1, 1].axvline(
        fixed["median_coherent_x_response"], color="#9b2c65", lw=2.5,
        label=f"event #12 median={fixed['median_coherent_x_response']:.3f}",
    )
    axes[1, 1].set(
        title=(f"#12 fixed p={comparison['event12']['fixed_p_raw']:.3f}; "
               f"joint p={comparison['event12']['joint_p_raw']:.3f}; "
               f"pair-shift p={comparison['fixed_pair_shift_p_raw']:.3f}"),
        xlabel="Median coherent X response", ylabel="Unrelated anchors",
    )
    axes[1, 1].legend()
    axes[1, 1].grid(axis="y", alpha=0.2)
    fig.suptitle("Two-candidate frozen comparison: event #12 follow-up", fontsize=15)
    path = OUT / "pds_20080113_2137_event12_2d_xfront_two_candidate.png"
    fig.savefig(path, dpi=230)
    plt.close(fig)
    return path


def write_report(comparison):
    e9, e12 = comparison["event9"], comparison["event12"]
    lines = [
        "# Two-candidate 2-D X-front comparison: events #9 and #12",
        "",
        "## Direct answer",
        "",
        f"If the event family were genuinely fixed to two independent candidates, the smaller "
        f"raw joint p={comparison['minimum_joint_p_raw']:.4f} would give a two-candidate "
        f"probability 1-(1-p)^2={comparison['two_candidate_probability']:.4f}.",
        "",
        "This does not retroactively replace the 12-event discovery scan, because event #9 was "
        "recognized as the strongest axial morphology after all 12 strong COR2 anchors were "
        "examined.  The m=2 number is therefore a descriptive candidate-family calculation.  "
        "It becomes confirmatory only on a new data representation or independent sample.",
        "",
        "## Results (raw, no BH)",
        "",
        "| Event | Frozen nodes | Positive coherent nodes | Fixed-template p | Joint p | Pass |",
        "| ---: | ---: | ---: | ---: | ---: | --- |",
        f"| #9 | {e9['node_count']} | {e9['positive_node_count']} | {e9['fixed_p_raw']:.4f} | "
        f"{e9['joint_p_raw']:.4f} | {'yes' if e9['passed'] else 'no'} |",
        f"| #12 | {e12['node_count']} | {e12['positive_node_count']} | {e12['fixed_p_raw']:.4f} | "
        f"{e12['joint_p_raw']:.4f} | {'yes' if e12['passed'] else 'no'} |",
        "",
        "An event is not rescued by a small joint p value if its fixed-template rarity or "
        "multi-node coherence condition fails.",
        "",
        "## Fixed 360-min pair-shift test",
        "",
        f"The real #9 -> #12 separation is exactly {comparison['pair_separation_min']:.0f} min.  "
        f"Both thresholds were applied after shifting the pair together through "
        f"{comparison['valid_shifted_pair_count']} valid control positions.  "
        f"{comparison['shifted_pair_exceedance_count']} shifted pairs equalled or exceeded both "
        f"events, giving raw pair p={comparison['fixed_pair_shift_p_raw']:.4f}.",
        "",
        "This is stronger than the minimum-p two-candidate calculation because it tests the "
        "ordered physical pair and preserves its separation.  It is nevertheless exploratory: "
        "the two candidates and this pair statistic were defined within the already inspected "
        "11--14 January interval, and the calibration is still Level-0-derived.",
        "",
        "## Event #12 node responses",
        "",
        "| Node | Vertex | |slope| | tB X | pB X | coherent X |",
        "| ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in comparison["event12_node_responses"]:
        lines.append(
            f"| {row['nominal_radius_rsun']:.3f} | {row['vertex_radius_rsun']:.3f} | "
            f"{row['absolute_slope_rsun_per_degree']:.4f} | {row['tb_x_response']:.3f} | "
            f"{row['pb_x_response']:.3f} | {row['coherent_x_response']:.3f} |"
        )
    lines += [
        "",
        "## Timing and displacement",
        "",
        "No additional event-dependent time shift was fitted.  Each radius uses the exact COR2 "
        "anchor plus the frozen propagation delay.  The PA axis is fixed.  Only the preregistered "
        "+/-0.075 R_sun vertex tolerance and frozen slope grid are scanned.",
        "",
        "## Limitation",
        "",
        "This remains a bias/exposure-normalized Level-0 fitpol morphology test.  The next truly "
        "confirmatory m=2 step is to apply the frozen pair to SECCHI_PREP-calibrated products and "
        "independent background representations.",
    ]
    path = OUT / "PDS_20080113_two_candidate_2D_Xfront_no_BH.md"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def main():
    grid, offset_cubes, maps = prepare_offset_cubes()
    radii = maps["radii"]
    raw_zmaps = {name: cells.robust_standardize(raw) for name, raw in maps["raw"].items()}
    families = core.model_families(radii)
    models = []
    for family in ("constant", "acceleration", "deceleration"):
        models.extend(families[family])
    tau_matrix = np.stack([model.tau for model in models])
    phase_results = json.loads(x2d.PHASE_RESULTS.read_text())
    reference = phase.model_from_dict(phase_results["reference_model"], radii)
    margin = abs(phase.PHASE_GRID[0]) + abs(phase.LOCAL_GRID[0])
    lo = maps["grid"][0] - np.min(reference.tau) + margin
    hi = maps["grid"][-1] - margin
    field_anchors = np.arange(
        math.ceil(lo / 30.0) * 30.0,
        math.floor(hi / 30.0) * 30.0 + 0.1,
        30.0,
    )

    anchor = float(maps["event_times"][CANDIDATE12["event_index"]])
    fixed_model = phase.model_from_dict(CANDIDATE12["model"], radii)
    event_maps = {
        product: x2d.aligned_slice(cube, grid, anchor, fixed_model.tau)
        for product, cube in offset_cubes.items()
    }
    event_fixed = score_top_n(
        event_maps, radii, CANDIDATE12["nodes_rsun"], len(CANDIDATE12["nodes_rsun"])
    )
    event_pipeline_info = cells.analyze_anchor(
        anchor, maps["grid"], radii, maps["zmaps"], raw_zmaps, models, tau_matrix
    )
    event_pipeline_model = phase.model_from_dict(event_pipeline_info["best_model"], radii)
    pipeline_maps = {
        product: x2d.aligned_slice(cube, grid, anchor, event_pipeline_model.tau)
        for product, cube in offset_cubes.items()
    }
    event_pipeline_nodes = np.asarray(
        event_pipeline_info["configs"]["strict"]["matched_primary_radii_rsun"], float
    )
    event_pipeline = score_top_n(
        pipeline_maps, radii, event_pipeline_nodes, len(CANDIDATE12["nodes_rsun"])
    )
    event_cell_stat = float(event_pipeline_info["configs"]["strict"]["cell_stat"])

    null_rows = []
    for number, null_anchor in enumerate(field_anchors, 1):
        fixed_maps = {
            product: x2d.aligned_slice(cube, grid, float(null_anchor), fixed_model.tau)
            for product, cube in offset_cubes.items()
        }
        fixed_score = score_top_n(
            fixed_maps, radii, CANDIDATE12["nodes_rsun"], len(CANDIDATE12["nodes_rsun"])
        )
        info = cells.analyze_anchor(
            float(null_anchor), maps["grid"], radii, maps["zmaps"], raw_zmaps, models, tau_matrix
        )
        model = phase.model_from_dict(info["best_model"], radii)
        selected_maps = {
            product: x2d.aligned_slice(cube, grid, float(null_anchor), model.tau)
            for product, cube in offset_cubes.items()
        }
        nodes = np.asarray(info["configs"]["strict"]["matched_primary_radii_rsun"], float)
        pipeline_score = score_top_n(
            selected_maps, radii, nodes, len(CANDIDATE12["nodes_rsun"])
        )
        null_rows.append({
            "anchor_min": float(null_anchor),
            "utc": anchor_utc(null_anchor),
            "fixed_median_x": fixed_score["median_coherent_x_response"],
            "fixed_positive_nodes": fixed_score["positive_node_count"],
            "pipeline_cell_stat": float(info["configs"]["strict"]["cell_stat"]),
            "pipeline_matched_peak_count": int(info["configs"]["strict"]["matched_peak_count"]),
            "pipeline_median_top3_x": pipeline_score["median_coherent_x_response"],
            "best_product": info["best_product"],
            "best_family": info["best_model"]["family"],
        })
        if number % 25 == 0:
            print(f"Event #12 2-D null: {number}/{len(field_anchors)}", flush=True)

    fixed_null = np.array([row["fixed_median_x"] for row in null_rows], float)
    fixed_p = empirical_p(event_fixed["median_coherent_x_response"], fixed_null)
    joint_flags = np.array([
        row["pipeline_cell_stat"] >= event_cell_stat
        and np.isfinite(row["pipeline_median_top3_x"])
        and row["pipeline_median_top3_x"] >= event_pipeline["median_coherent_x_response"]
        for row in null_rows
    ], bool)
    joint_k = int(np.sum(joint_flags))
    joint_p = float((joint_k + 1) / (len(null_rows) + 1))
    for row, flag in zip(null_rows, joint_flags):
        row["joint_exceedance"] = bool(flag)

    e12_pass = bool(
        fixed_p < 0.05 and joint_p < 0.05
        and event_fixed["positive_node_count"] >= CANDIDATE12["required_positive_nodes"]
    )
    e9_full = json.loads(EVENT9_RESULTS.read_text())
    e9_fixed = e9_full["event9_fixed_template"]
    e9_pass = bool(e9_full["passed_preregistered_detection_rule"])
    e9_joint = float(e9_full["joint_p_raw"])
    min_joint = min(e9_joint, joint_p)
    p_two = float(1.0 - (1.0 - min_joint) ** 2)

    with EVENT9_NULL_TABLE.open(newline="", encoding="utf-8") as handle:
        event9_null = {float(row["anchor_min"]): row for row in csv.DictReader(handle)}
    event12_null = {float(row["anchor_min"]): row for row in null_rows}
    event9_cell_threshold = float(e9_full["event9_pipeline"]["cell_stat"])
    event9_x_threshold = float(
        e9_full["event9_pipeline"]["median_coherent_x_response"]
    )

    def finite_ge(value, threshold):
        try:
            numeric = float(value)
        except (TypeError, ValueError):
            return False
        return bool(np.isfinite(numeric) and numeric >= threshold)

    pair_rows = []
    for early_anchor, early in sorted(event9_null.items()):
        late_anchor = early_anchor + PAIR_SEPARATION_MIN
        late = event12_null.get(late_anchor)
        if late is None:
            continue
        early_pass = (
            finite_ge(early["pipeline_cell_stat"], event9_cell_threshold)
            and finite_ge(early["pipeline_median_top4_x"], event9_x_threshold)
        )
        late_pass = (
            finite_ge(late["pipeline_cell_stat"], event_cell_stat)
            and finite_ge(
                late["pipeline_median_top3_x"],
                event_pipeline["median_coherent_x_response"],
            )
        )
        pair_rows.append({
            "early_anchor_min": early_anchor,
            "early_utc": early["utc"],
            "late_anchor_min": late_anchor,
            "late_utc": late["utc"],
            "event9_like_pass": early_pass,
            "event12_like_pass": late_pass,
            "pair_pass": bool(early_pass and late_pass),
            "early_cell_stat": early["pipeline_cell_stat"],
            "early_median_top4_x": early["pipeline_median_top4_x"],
            "late_cell_stat": late["pipeline_cell_stat"],
            "late_median_top3_x": late["pipeline_median_top3_x"],
        })
    pair_k = int(sum(row["pair_pass"] for row in pair_rows))
    pair_p = float((pair_k + 1) / (len(pair_rows) + 1))

    comparison = {
        "candidate_family_size": 2,
        "two_candidate_status": "descriptive; event #9 was selected after the 12-anchor scan",
        "event9": {
            "node_count": 4,
            "positive_node_count": int(e9_fixed["positive_node_count"]),
            "fixed_p_raw": float(e9_fixed["p_raw"]),
            "joint_p_raw": e9_joint,
            "passed": e9_pass,
        },
        "event12": {
            "node_count": 3,
            "positive_node_count": int(event_fixed["positive_node_count"]),
            "fixed_median_coherent_x_response": float(event_fixed["median_coherent_x_response"]),
            "pipeline_median_coherent_x_response": float(event_pipeline["median_coherent_x_response"]),
            "cell_stat": event_cell_stat,
            "fixed_p_raw": fixed_p,
            "joint_null_exceedances": joint_k,
            "joint_p_raw": joint_p,
            "passed": e12_pass,
        },
        "event12_node_responses": event_fixed["nodes"],
        "minimum_joint_p_raw": min_joint,
        "two_candidate_probability": p_two,
        "two_candidate_family_below_0_05": bool(p_two < 0.05),
        "pair_separation_min": PAIR_SEPARATION_MIN,
        "valid_shifted_pair_count": int(len(pair_rows)),
        "shifted_pair_exceedance_count": pair_k,
        "fixed_pair_shift_p_raw": pair_p,
        "physical_detection": bool(p_two < 0.05 and (e9_pass or e12_pass)),
        "null_anchor_count": int(len(null_rows)),
        "no_BH": True,
        "no_additional_time_shift": True,
    }

    figure_path = make_figure(event_maps, radii, event_fixed, fixed_null, comparison)
    report_path = write_report(comparison)
    comparison["files"] = {"figure": figure_path.name, "report": report_path.name}
    json_path = OUT / "pds_20080113_two_candidate_2d_xfront_results_no_BH.json"
    json_path.write_text(json.dumps(comparison, indent=2), encoding="utf-8")
    csv_path = OUT / "pds_20080113_event12_2d_xfront_null_table.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(null_rows[0]))
        writer.writeheader()
        writer.writerows(null_rows)
    pair_csv_path = OUT / "pds_20080113_event09_event12_fixed_360min_pair_null.csv"
    with pair_csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(pair_rows[0]))
        writer.writeheader()
        writer.writerows(pair_rows)
    np.savez_compressed(
        OUT / "pds_20080113_2137_event12_2d_xfront_maps.npz",
        radii=radii, pa_offsets=x2d.PA_OFFSETS,
        event_total_b=event_maps["total_b"], event_pb=event_maps["pb"],
        fixed_null_median_x=fixed_null,
    )
    print(json.dumps(comparison, indent=2))


if __name__ == "__main__":
    main()
