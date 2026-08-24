#!/usr/bin/env python3
"""Frozen time-resolved X-node activation-order test for events #9 and #12.

The specification is stored in
PDS_20080113_node_activation_order_preregistered.md and was written before
the activation curves were computed.
"""

from __future__ import annotations

import csv
import json
import math
from datetime import timedelta
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import analyze_event09_2d_xfront_20080113 as x2d
import analyze_pds_event_phase_jitter_20080111_14 as phase
import analyze_pds_nonlinear_transport_20080111_14 as core
import analyze_two_candidate_2d_xfront_20080113 as pair_test
import analyze_two_event_representation_robustness_20080113 as reps


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "results" / "pds_20080111_14_nonlinear_level0"
OUT = SOURCE / "node_activation_order"
OUT.mkdir(parents=True, exist_ok=True)

PAIR_NULL_FILE = (
    SOURCE / "two_candidate_2d_xfront"
    / "pds_20080113_event09_event12_fixed_360min_pair_null.csv"
)

REPRESENTATIONS = ("bff", "base60", "nrgf60")
OFFSETS_MIN = np.arange(-45.0, 45.1, core.CADENCE_MIN)
ZERO_INDEX = int(np.argmin(np.abs(OFFSETS_MIN)))
RSUN_KM = 695700.0
PAIR_SEPARATION_MIN = 360.0

CANDIDATES = {
    "event9": {
        "event_number": 9,
        "event_index": 8,
        "utc": core.STRONG_EVENTS[8][0],
        "nodes": np.array([1.775, 2.125, 2.500, 2.825], float),
        "model": {
            "family": "deceleration",
            "v_inner_km_s": 100.0,
            "v_outer_km_s": 30.0,
            "change_radius_rsun": 2.20,
        },
    },
    "event12": {
        "event_number": 12,
        "event_index": 11,
        "utc": core.STRONG_EVENTS[11][0],
        "nodes": np.array([1.950, 2.350, 2.700], float),
        "model": {
            "family": "deceleration",
            "v_inner_km_s": 100.0,
            "v_outer_km_s": 25.0,
            "change_radius_rsun": 2.90,
        },
    },
}


def empirical_p(exceedances, total):
    return float((int(exceedances) + 1) / (int(total) + 1))


def anchor_utc(anchor):
    return (core.T0 + timedelta(minutes=float(anchor))).strftime("%Y-%m-%dT%H:%M:%SZ")


def fixed_line_response(map2d, radii, geometry):
    vertex = float(geometry["vertex_radius_rsun"])
    slope = float(geometry["absolute_slope_rsun_per_degree"])
    positive = x2d.diagonal_response(map2d, radii, vertex, slope)
    negative = x2d.diagonal_response(map2d, radii, vertex, -slope)
    if not (np.isfinite(positive) and np.isfinite(negative)):
        return np.nan
    return float(min(positive, negative))


def robust_zero_prominence(curve):
    curve = np.asarray(curve, float)
    valid = np.isfinite(curve)
    if not valid[ZERO_INDEX] or np.sum(valid) < 4:
        return np.nan
    center = float(np.nanmedian(curve))
    mad = float(1.4826 * np.nanmedian(np.abs(curve[valid] - center)))
    std = float(np.nanstd(curve[valid]))
    scale = mad if mad > 1e-8 else (std if std > 1e-8 else 1.0)
    return float((curve[ZERO_INDEX] - center) / scale)


def peak_offset(curve):
    curve = np.asarray(curve, float)
    if not np.any(np.isfinite(curve)):
        return np.nan, np.nan
    index = int(np.nanargmax(curve))
    return float(OFFSETS_MIN[index]), float(curve[index])


def activation_result(offset_products, grid, anchor, model, radii, nodes, geometries):
    curves = {
        "total_b": np.full((len(nodes), len(OFFSETS_MIN)), np.nan),
        "pb": np.full((len(nodes), len(OFFSETS_MIN)), np.nan),
    }
    for k, offset in enumerate(OFFSETS_MIN):
        maps = {
            product: x2d.aligned_slice(cube, grid, float(anchor + offset), model.tau)
            for product, cube in offset_products.items()
        }
        for i, geometry in enumerate(geometries):
            for product in ("total_b", "pb"):
                curves[product][i, k] = fixed_line_response(
                    maps[product], radii, geometry
                )

    coherent = np.minimum(curves["total_b"], curves["pb"])
    tb_offsets, pb_offsets, coherent_offsets = [], [], []
    tb_peaks, pb_peaks, coherent_peaks = [], [], []
    node_prominence = []
    for i in range(len(nodes)):
        offset, peak = peak_offset(curves["total_b"][i])
        tb_offsets.append(offset)
        tb_peaks.append(peak)
        offset, peak = peak_offset(curves["pb"][i])
        pb_offsets.append(offset)
        pb_peaks.append(peak)
        offset, peak = peak_offset(coherent[i])
        coherent_offsets.append(offset)
        coherent_peaks.append(peak)
        node_prominence.append(robust_zero_prominence(coherent[i]))

    tb_offsets = np.asarray(tb_offsets, float)
    pb_offsets = np.asarray(pb_offsets, float)
    coherent_offsets = np.asarray(coherent_offsets, float)
    tau_nodes = np.interp(nodes, radii, model.tau)
    activation_times = float(anchor) + tau_nodes + coherent_offsets

    segment_speeds = []
    segment_labels = []
    for i in range(len(nodes) - 1):
        dt_min = activation_times[i + 1] - activation_times[i]
        speed = (nodes[i + 1] - nodes[i]) * RSUN_KM / (dt_min * 60.0)
        segment_speeds.append(float(speed))
        segment_labels.append(f"{nodes[i]:.3f}-{nodes[i+1]:.3f}")
    dt_min = float(anchor) - activation_times[-1]
    speed = (3.0 - nodes[-1]) * RSUN_KM / (dt_min * 60.0)
    segment_speeds.append(float(speed))
    segment_labels.append(f"{nodes[-1]:.3f}-3.000")
    segment_speeds = np.asarray(segment_speeds, float)

    required = int(math.ceil(0.75 * len(nodes)))
    positive_zero = int(np.sum(coherent[:, ZERO_INDEX] > 0))
    agreement = int(np.sum(np.abs(tb_offsets - pb_offsets) <= core.CADENCE_MIN))
    ordered = bool(np.all(np.diff(activation_times) > 0))
    speeds_valid = bool(np.all(
        np.isfinite(segment_speeds)
        & (segment_speeds >= 25.0)
        & (segment_speeds <= 300.0)
    ))
    boundary_count = int(np.sum(np.abs(coherent_offsets) == np.max(np.abs(OFFSETS_MIN))))
    boundary_allowed = int(math.floor(0.25 * len(nodes)))
    quality = bool(
        positive_zero >= required
        and agreement >= required
        and ordered
        and speeds_valid
        and boundary_count <= boundary_allowed
    )

    node_rows = []
    for i, node in enumerate(nodes):
        node_rows.append({
            "radius_rsun": float(node),
            "vertex_radius_rsun": float(geometries[i]["vertex_radius_rsun"]),
            "absolute_slope_rsun_per_degree": float(
                geometries[i]["absolute_slope_rsun_per_degree"]
            ),
            "tb_peak_offset_min": float(tb_offsets[i]),
            "pb_peak_offset_min": float(pb_offsets[i]),
            "coherent_peak_offset_min": float(coherent_offsets[i]),
            "tb_peak_x": float(tb_peaks[i]),
            "pb_peak_x": float(pb_peaks[i]),
            "coherent_peak_x": float(coherent_peaks[i]),
            "coherent_x_at_zero": float(coherent[i, ZERO_INDEX]),
            "temporal_prominence_z": float(node_prominence[i]),
            "frozen_tau_min": float(tau_nodes[i]),
            "activation_min_from_t0": float(activation_times[i]),
            "activation_utc": anchor_utc(activation_times[i]),
        })

    return {
        "score": float(np.nanmedian(node_prominence)),
        "node_prominence_z": np.asarray(node_prominence, float),
        "curves": {**curves, "coherent": coherent},
        "nodes": node_rows,
        "positive_zero_nodes": positive_zero,
        "required_nodes": required,
        "tb_pb_agreement_nodes": agreement,
        "ordered_activation_times": ordered,
        "segment_labels": segment_labels,
        "segment_speeds_km_s": segment_speeds,
        "all_segment_speeds_valid": speeds_valid,
        "boundary_peak_nodes": boundary_count,
        "boundary_allowed": boundary_allowed,
        "quality_pass": quality,
    }


def prepare_inputs():
    bff_grid, bff_offsets, maps = pair_test.prepare_offset_cubes()
    grid, radii_full, radii, pa, path_full, event_times, regular = reps.regularize_cube()
    if not (np.array_equal(grid, bff_grid) and np.allclose(radii, maps["radii"])):
        raise RuntimeError("Primary and alternative-representation grids do not match")
    derived = reps.build_representations(regular, radii_full, radii, pa, path_full)
    products = {
        "bff": bff_offsets,
        "base60": derived["base60"],
        "nrgf60": derived["nrgf60"],
    }
    return grid, radii, event_times, products


def observed_geometries(products, grid, radii, event_times, models):
    geometries = {name: {} for name in CANDIDATES}
    zero_maps = {name: {} for name in CANDIDATES}
    for candidate_name, candidate in CANDIDATES.items():
        anchor = float(event_times[candidate["event_index"]])
        for representation in REPRESENTATIONS:
            maps = {
                product: x2d.aligned_slice(cube, grid, anchor, models[candidate_name].tau)
                for product, cube in products[representation].items()
            }
            score = pair_test.score_top_n(
                maps, radii, candidate["nodes"], len(candidate["nodes"])
            )
            geometries[candidate_name][representation] = score["nodes"]
            zero_maps[candidate_name][representation] = maps
    return geometries, zero_maps


def serializable_result(result):
    return {
        key: value for key, value in result.items()
        if key not in ("curves", "node_prominence_z", "segment_speeds_km_s")
    } | {
        "node_prominence_z": result["node_prominence_z"].tolist(),
        "segment_speeds_km_s": result["segment_speeds_km_s"].tolist(),
    }


def make_activation_figure(observed):
    fig, axes = plt.subplots(3, 2, figsize=(13.5, 13.5), constrained_layout=True)
    colors = ("#4c78a8", "#f58518", "#54a24b", "#b279a2")
    for row, representation in enumerate(REPRESENTATIONS):
        for col, candidate_name in enumerate(("event9", "event12")):
            candidate = CANDIDATES[candidate_name]
            result = observed[candidate_name][representation]
            ax = axes[row, col]
            for i, node in enumerate(candidate["nodes"]):
                ax.plot(
                    OFFSETS_MIN, result["curves"]["coherent"][i], marker="o",
                    color=colors[i], label=f"{node:.3f} R_sun",
                )
                peak = result["nodes"][i]["coherent_peak_offset_min"]
                ax.axvline(peak, color=colors[i], alpha=0.22, lw=1)
            ax.axvline(0, color="black", ls="--", lw=1.2, label="frozen path")
            ax.axhline(0, color="0.5", ls=":", lw=0.8)
            ax.set(
                title=(f"{representation} | event #{candidate['event_number']} | "
                       f"score={result['score']:.2f} | "
                       f"quality={'yes' if result['quality_pass'] else 'no'}"),
                xlabel="Offset from frozen propagation curve (min)",
                ylabel="Coherent min(tB, pB) X response",
            )
            ax.legend(fontsize=8, ncol=2)
            ax.grid(alpha=0.2)
    path = OUT / "pds_20080113_node_activation_curves.png"
    fig.savefig(path, dpi=230)
    plt.close(fig)
    return path


def make_timing_figure(observed, event_times, models, radii):
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.5), constrained_layout=True)
    colors = {"bff": "#4c78a8", "base60": "#f58518", "nrgf60": "#54a24b"}
    for ax, candidate_name in zip(axes, ("event9", "event12")):
        candidate = CANDIDATES[candidate_name]
        anchor = float(event_times[candidate["event_index"]])
        nodes = candidate["nodes"]
        expected = anchor + np.interp(nodes, radii, models[candidate_name].tau)
        ax.plot(expected - anchor, nodes, color="black", ls="--", lw=2, label="frozen model")
        ax.scatter([0], [3.0], marker="*", s=90, color="black", zorder=4, label="COR2 anchor")
        for representation in REPRESENTATIONS:
            times = np.array([
                node["activation_min_from_t0"]
                for node in observed[candidate_name][representation]["nodes"]
            ])
            ax.plot(times - anchor, nodes, marker="o", color=colors[representation],
                    label=representation)
        ax.set(
            title=f"Event #{candidate['event_number']} activation order",
            xlabel="Minutes relative to COR2 3-R_sun anchor",
            ylabel="Radius (R_sun)",
        )
        ax.legend(fontsize=8)
        ax.grid(alpha=0.2)
    path = OUT / "pds_20080113_node_activation_times.png"
    fig.savefig(path, dpi=230)
    plt.close(fig)
    return path


def write_report(summary):
    lines = [
        "# Time-resolved node-activation order: events #9 and #12",
        "",
        "## Verdict",
        "",
        summary["verdict"],
        "",
        "## Frozen event results (raw; no BH)",
        "",
        "| Event | Representation | Score | Positive at zero | tB/pB agreement | "
        "Ordered | Speeds valid | Boundary peaks | Quality | Formal score-null p | Eligible |",
        "| ---: | --- | ---: | ---: | ---: | --- | --- | ---: | --- | ---: | --- |",
    ]
    for candidate_name in ("event9", "event12"):
        event_number = CANDIDATES[candidate_name]["event_number"]
        for representation in REPRESENTATIONS:
            row = summary["events"][candidate_name][representation]
            lines.append(
                f"| #{event_number} | {representation} | {row['score']:.3f} | "
                f"{row['positive_zero_nodes']}/{row['node_count']} | "
                f"{row['tb_pb_agreement_nodes']}/{row['node_count']} | "
                f"{'yes' if row['ordered_activation_times'] else 'no'} | "
                f"{'yes' if row['all_segment_speeds_valid'] else 'no'} | "
                f"{row['boundary_peak_nodes']} | "
                f"{'pass' if row['quality_pass'] else 'fail'} | "
                f"{row['individual_formal_p_raw']:.4f} | "
                f"{'yes' if row['individual_p_eligible'] else 'no'} |"
            )
    lines += [
        "",
        "## Fixed 360-min pair null",
        "",
        "| Representation | Quality controls early/late/both | Pair exceedances | "
        "Valid shifted pairs | Formal raw p | Eligible |",
        "| --- | ---: | ---: | ---: | ---: | --- |",
    ]
    for representation in REPRESENTATIONS:
        row = summary["pair_null"][representation]
        lines.append(
            f"| {representation} | {row['early_quality_controls']}/"
            f"{row['late_quality_controls']}/{row['both_quality_controls']} | "
            f"{row['exceedances']} | {row['valid_pairs']} | {row['formal_p_raw']:.4f} | "
            f"{'yes' if row['p_eligible'] else 'no'} |"
        )
    lines += [
        "",
        f"Primary >=2-of-3 pair exceedances: {summary['primary_pair_exceedances']}/"
        f"{summary['valid_pair_count']}; formal raw p="
        f"{summary['primary_pair_formal_p_raw']:.4f}; eligible: "
        f"{'yes' if summary['primary_pair_p_eligible'] else 'no'}.",
        "",
        f"Event #9 passes the frozen quality rule in "
        f"{summary['event9_quality_representation_count']}/3 representations; "
        f"event #12 in {summary['event12_quality_representation_count']}/3.",
        "",
        "The formal 1/(N+1) values above are not detection probabilities when the real "
        "event fails the frozen quality rule. Here both targets are ineligible, so the "
        "0.0076 floor cannot rescue the activation-order claim.",
        "",
        "## Node timings and segment speeds",
        "",
    ]
    for candidate_name in ("event9", "event12"):
        event_number = CANDIDATES[candidate_name]["event_number"]
        for representation in REPRESENTATIONS:
            row = summary["events"][candidate_name][representation]
            lines += [
                f"### Event #{event_number}, {representation}",
                "",
                "| R_sun | tB peak (min) | pB peak (min) | Coherent peak (min) | "
                "X at zero | Activation UTC |",
                "| ---: | ---: | ---: | ---: | ---: | --- |",
            ]
            for node in row["nodes"]:
                lines.append(
                    f"| {node['radius_rsun']:.3f} | {node['tb_peak_offset_min']:.0f} | "
                    f"{node['pb_peak_offset_min']:.0f} | "
                    f"{node['coherent_peak_offset_min']:.0f} | "
                    f"{node['coherent_x_at_zero']:.3f} | {node['activation_utc']} |"
                )
            speeds = ", ".join(
                f"{label}: {speed:.1f} km/s"
                for label, speed in zip(row["segment_labels"], row["segment_speeds_km_s"])
            )
            lines += ["", f"Segment pattern speeds: {speeds}.", ""]
    lines += [
        "## Interpretation boundary",
        "",
        "The test distinguishes time-ordered propagation on the frozen path from a simultaneous "
        "or incoherent X morphology. It does not establish a local density compression, a "
        "Rankine--Hugoniot solution, a magnetosonic Mach number, or a slow/fast MHD branch. "
        "The inputs remain Level-0-derived COR1 fitpol diagnostics, and formal SolarSoft "
        "SECCHI_PREP Level-1 confirmation is still required for a shock claim.",
        "",
        "All probabilities are raw; no BH correction is used.",
    ]
    path = OUT / "PDS_20080113_node_activation_order_no_BH.md"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def main():
    grid, radii, event_times, products = prepare_inputs()
    models = {
        name: phase.model_from_dict(candidate["model"], radii)
        for name, candidate in CANDIDATES.items()
    }
    geometries, _ = observed_geometries(products, grid, radii, event_times, models)

    observed = {name: {} for name in CANDIDATES}
    for candidate_name, candidate in CANDIDATES.items():
        anchor = float(event_times[candidate["event_index"]])
        for representation in REPRESENTATIONS:
            observed[candidate_name][representation] = activation_result(
                products[representation], grid, anchor, models[candidate_name], radii,
                candidate["nodes"], geometries[candidate_name][representation],
            )

    with PAIR_NULL_FILE.open(newline="", encoding="utf-8") as handle:
        early_anchors = np.array([
            float(row["early_anchor_min"]) for row in csv.DictReader(handle)
        ])

    null_rows = []
    for number, early_anchor in enumerate(early_anchors, 1):
        late_anchor = float(early_anchor + PAIR_SEPARATION_MIN)
        row = {
            "early_anchor_min": float(early_anchor),
            "late_anchor_min": late_anchor,
            "early_utc": anchor_utc(early_anchor),
            "late_utc": anchor_utc(late_anchor),
        }
        rep_exceedance_count = 0
        fully_valid = True
        for representation in REPRESENTATIONS:
            early = activation_result(
                products[representation], grid, early_anchor, models["event9"], radii,
                CANDIDATES["event9"]["nodes"], geometries["event9"][representation],
            )
            late = activation_result(
                products[representation], grid, late_anchor, models["event12"], radii,
                CANDIDATES["event12"]["nodes"], geometries["event12"][representation],
            )
            valid = bool(np.isfinite(early["score"]) and np.isfinite(late["score"]))
            early_exceed = bool(
                valid and early["quality_pass"]
                and early["score"] >= observed["event9"][representation]["score"]
            )
            late_exceed = bool(
                valid and late["quality_pass"]
                and late["score"] >= observed["event12"][representation]["score"]
            )
            pair_exceed = bool(early_exceed and late_exceed)
            row[f"{representation}_early_score"] = float(early["score"])
            row[f"{representation}_late_score"] = float(late["score"])
            row[f"{representation}_early_quality"] = bool(early["quality_pass"])
            row[f"{representation}_late_quality"] = bool(late["quality_pass"])
            row[f"{representation}_early_exceed"] = early_exceed
            row[f"{representation}_late_exceed"] = late_exceed
            row[f"{representation}_pair_exceed"] = pair_exceed
            fully_valid &= valid
            rep_exceedance_count += int(pair_exceed)
        row["fully_valid"] = fully_valid
        row["representation_pair_exceedance_count"] = rep_exceedance_count
        row["primary_pair_exceed"] = bool(fully_valid and rep_exceedance_count >= 2)
        null_rows.append(row)
        if number % 20 == 0:
            print(f"Activation-order pair null: {number}/{len(early_anchors)}", flush=True)

    valid_rows = [row for row in null_rows if row["fully_valid"]]
    pair_null = {}
    event_p = {name: {} for name in CANDIDATES}
    for representation in REPRESENTATIONS:
        pair_k = int(sum(row[f"{representation}_pair_exceed"] for row in valid_rows))
        early_k = int(sum(row[f"{representation}_early_exceed"] for row in valid_rows))
        late_k = int(sum(row[f"{representation}_late_exceed"] for row in valid_rows))
        early_quality_controls = int(sum(
            row[f"{representation}_early_quality"] for row in valid_rows
        ))
        late_quality_controls = int(sum(
            row[f"{representation}_late_quality"] for row in valid_rows
        ))
        both_quality_controls = int(sum(
            row[f"{representation}_early_quality"]
            and row[f"{representation}_late_quality"] for row in valid_rows
        ))
        p_eligible = bool(
            observed["event9"][representation]["quality_pass"]
            and observed["event12"][representation]["quality_pass"]
        )
        pair_null[representation] = {
            "exceedances": pair_k,
            "valid_pairs": len(valid_rows),
            "formal_p_raw": empirical_p(pair_k, len(valid_rows)),
            "p_eligible": p_eligible,
            "early_quality_controls": early_quality_controls,
            "late_quality_controls": late_quality_controls,
            "both_quality_controls": both_quality_controls,
        }
        event_p["event9"][representation] = empirical_p(early_k, len(valid_rows))
        event_p["event12"][representation] = empirical_p(late_k, len(valid_rows))

    primary_k = int(sum(row["primary_pair_exceed"] for row in valid_rows))
    primary_p = empirical_p(primary_k, len(valid_rows))
    quality_counts = {
        name: int(sum(observed[name][rep]["quality_pass"] for rep in REPRESENTATIONS))
        for name in CANDIDATES
    }
    passed = bool(
        primary_p < 0.05
        and quality_counts["event9"] >= 2
        and quality_counts["event12"] >= 2
    )
    primary_p_eligible = bool(
        quality_counts["event9"] >= 2 and quality_counts["event12"] >= 2
    )
    if passed:
        verdict = (
            "The ordered #9 -> #12 pair passes the frozen time-resolved activation-order "
            "rule. Both events retain tB/pB-consistent outward node ordering in at least two "
            "representations, and the >=2-of-3 shifted-pair null is locally unusual. This "
            "supports pattern transport through the X/diamond-like morphology, but not a "
            "physical shock identification."
        )
    else:
        verdict = (
            "The ordered #9 -> #12 pair does not pass the complete frozen time-resolved "
            "activation-order rule. The X/diamond-like morphology remains processing-robust, "
            "but the present data do not show a reproducible tB/pB-consistent outward "
            "activation sequence through the frozen nodes."
        )

    summary_events = {name: {} for name in CANDIDATES}
    for candidate_name, candidate in CANDIDATES.items():
        for representation in REPRESENTATIONS:
            result = serializable_result(observed[candidate_name][representation])
            result["node_count"] = len(candidate["nodes"])
            result["individual_formal_p_raw"] = event_p[candidate_name][representation]
            result["individual_p_eligible"] = bool(result["quality_pass"])
            summary_events[candidate_name][representation] = result

    summary = {
        "events": summary_events,
        "fixed_candidates": {
            name: {
                "event_number": candidate["event_number"],
                "utc": candidate["utc"],
                "nodes_rsun": candidate["nodes"].tolist(),
                "model": candidate["model"],
            }
            for name, candidate in CANDIDATES.items()
        },
        "representations": list(REPRESENTATIONS),
        "offset_grid_min": OFFSETS_MIN.tolist(),
        "pair_separation_min": PAIR_SEPARATION_MIN,
        "pair_null": pair_null,
        "valid_pair_count": len(valid_rows),
        "primary_pair_exceedances": primary_k,
        "primary_pair_formal_p_raw": primary_p,
        "primary_pair_p_eligible": primary_p_eligible,
        "event9_quality_representation_count": quality_counts["event9"],
        "event12_quality_representation_count": quality_counts["event12"],
        "passed_frozen_activation_order_rule": passed,
        "verdict": verdict,
        "no_BH": True,
        "calibration": "bias/exposure-normalized Level-0-derived COR1 fitpol diagnostic",
    }

    curves_figure = make_activation_figure(observed)
    timing_figure = make_timing_figure(observed, event_times, models, radii)
    report_path = write_report(summary)
    summary["files"] = {
        "activation_curves": curves_figure.name,
        "activation_times": timing_figure.name,
        "report": report_path.name,
    }
    json_path = OUT / "pds_20080113_node_activation_order_results_no_BH.json"
    json_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    csv_path = OUT / "pds_20080113_node_activation_order_pair_null.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(null_rows[0]))
        writer.writeheader()
        writer.writerows(null_rows)

    payload = {"offsets_min": OFFSETS_MIN}
    for candidate_name in CANDIDATES:
        for representation in REPRESENTATIONS:
            result = observed[candidate_name][representation]
            for product in ("total_b", "pb", "coherent"):
                payload[f"{candidate_name}_{representation}_{product}"] = result["curves"][product]
    np.savez_compressed(OUT / "pds_20080113_node_activation_curves.npz", **payload)

    print(json.dumps({
        "passed": passed,
        "primary_pair_formal_p_raw": primary_p,
        "primary_pair_p_eligible": primary_p_eligible,
        "quality_counts": quality_counts,
        "pair_null": pair_null,
        "event_results": {
            name: {
                rep: {
                    "score": observed[name][rep]["score"],
                    "quality": observed[name][rep]["quality_pass"],
                    "positive_zero": observed[name][rep]["positive_zero_nodes"],
                    "agreement": observed[name][rep]["tb_pb_agreement_nodes"],
                    "ordered": observed[name][rep]["ordered_activation_times"],
                    "speeds": observed[name][rep]["segment_speeds_km_s"].tolist(),
                    "offsets": [n["coherent_peak_offset_min"] for n in observed[name][rep]["nodes"]],
                }
                for rep in REPRESENTATIONS
            }
            for name in CANDIDATES
        },
    }, indent=2))


if __name__ == "__main__":
    main()
