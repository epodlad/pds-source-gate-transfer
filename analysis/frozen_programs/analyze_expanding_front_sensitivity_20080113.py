#!/usr/bin/env python3
"""Physics-aware expanding-front sensitivity for COR1-A events #9 and #12."""

from __future__ import annotations

import csv
import json
import math
from datetime import timedelta
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import analyze_dynamic_geometry_gate as geometry
import analyze_node_activation_order_20080113 as strict
import analyze_pds_event_phase_jitter_20080111_14 as phase
import analyze_pds_nonlinear_transport_20080111_14 as core
import analyze_event09_2d_xfront_20080113 as x2d


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "results" / "pds_20080111_14_nonlinear_level0"
OUT = SOURCE / "expanding_front_sensitivity"
OUT.mkdir(parents=True, exist_ok=True)

OFFSETS_MIN = np.arange(-60.0, 60.1, core.CADENCE_MIN)
ANCHOR_JITTER_MIN = np.array([-15.0, 0.0, 15.0])
REPRESENTATIONS = strict.REPRESENTATIONS
CANDIDATES = strict.CANDIDATES
PAIR_NULL_FILE = strict.PAIR_NULL_FILE
RSUN_KM = strict.RSUN_KM
SPEED_RANGE = (25.0, 300.0)


def empirical_p(exceedances, total):
    return float((int(exceedances) + 1) / (int(total) + 1))


def utc_from_minute(value):
    return (core.T0 + timedelta(minutes=float(value))).strftime("%Y-%m-%dT%H:%M:%SZ")


def robust_probability(curve):
    curve = np.asarray(curve, float)
    valid = np.isfinite(curve)
    result = np.zeros_like(curve)
    if np.sum(valid) < 4:
        result[valid] = 1.0 / max(np.sum(valid), 1)
        return result
    center = float(np.nanmedian(curve))
    mad = float(1.4826 * np.nanmedian(np.abs(curve[valid] - center)))
    std = float(np.nanstd(curve[valid]))
    scale = mad if mad > 1e-8 else (std if std > 1e-8 else 1.0)
    z = np.clip((curve[valid] - center) / scale, -4.0, 4.0)
    weight = np.exp(z - np.max(z))
    result[valid] = weight / np.sum(weight)
    return result


def broaden_probability(probability, sigma_min):
    probability = np.asarray(probability, float)
    distance = OFFSETS_MIN[:, None] - OFFSETS_MIN[None, :]
    kernel = np.exp(-0.5 * np.square(distance / float(sigma_min)))
    kernel /= np.sum(kernel, axis=0, keepdims=True)
    result = kernel @ probability
    return result / np.sum(result)


def posterior_quantiles(probability):
    probability = np.asarray(probability, float)
    probability = probability / np.sum(probability)
    cdf = np.cumsum(probability)
    def quantile(q):
        return float(OFFSETS_MIN[min(int(np.searchsorted(cdf, q)), len(cdf) - 1)])
    mean = float(np.sum(probability * OFFSETS_MIN))
    mode = float(OFFSETS_MIN[int(np.argmax(probability))])
    variance = float(np.sum(probability * np.square(OFFSETS_MIN - mean)))
    return {
        "mode_min": mode,
        "mean_min": mean,
        "sigma_min": math.sqrt(max(variance, 0.0)),
        "q16_min": quantile(0.16),
        "q84_min": quantile(0.84),
    }


def activation_curves(offset_products, grid, anchor, model, radii, geometries):
    curves = {
        "total_b": np.full((len(geometries), len(OFFSETS_MIN)), np.nan),
        "pb": np.full((len(geometries), len(OFFSETS_MIN)), np.nan),
    }
    for k, offset in enumerate(OFFSETS_MIN):
        maps = {
            product: x2d.aligned_slice(cube, grid, float(anchor + offset), model.tau)
            for product, cube in offset_products.items()
        }
        for i, line in enumerate(geometries):
            for product in ("total_b", "pb"):
                curves[product][i, k] = strict.fixed_line_response(
                    maps[product], radii, line
                )
    return curves


def event_area(moment, anchor, model, radii, nodes):
    tau_full = np.interp(moment["radii"], radii, model.tau)
    targets = float(anchor) + tau_full
    log_area_path = geometry.interpolate_path(
        moment["log_area"], moment["times"], targets
    )
    good = np.isfinite(log_area_path)
    if np.sum(good) >= 3:
        node_log_area = np.interp(nodes, moment["radii"][good], log_area_path[good])
        measured = np.isfinite(node_log_area)
    else:
        node_log_area = np.full(len(nodes), np.nan)
        measured = np.zeros(len(nodes), bool)

    spherical_ratio = np.square(nodes / nodes[0])
    area_ratio = spherical_ratio.copy()
    if np.isfinite(node_log_area[0]):
        candidate_ratio = np.exp(node_log_area - node_log_area[0])
        use = np.isfinite(candidate_ratio) & (candidate_ratio > 0)
        area_ratio[use] = candidate_ratio[use]
        measured = use
    area_ratio[0] = 1.0
    sigma = np.clip(15.0 * np.sqrt(area_ratio), 15.0, 30.0)
    return {
        "area_ratio": area_ratio,
        "timing_sigma_min": sigma,
        "measured": measured,
        "node_log_area": node_log_area,
    }


def enumerate_sequence_probability(posteriors, anchor, tau_nodes, nodes):
    n_nodes, n_offsets = posteriors.shape
    meshes = np.meshgrid(*([np.arange(n_offsets)] * n_nodes), indexing="ij")
    indices = np.stack([mesh.ravel() for mesh in meshes], axis=1)
    offsets = OFFSETS_MIN[indices]
    weights = np.ones(len(indices), float)
    for i in range(n_nodes):
        weights *= posteriors[i, indices[:, i]]
    weights /= np.sum(weights)

    activation = float(anchor) + tau_nodes[None, :] + offsets
    ordered = np.all(np.diff(activation, axis=1) > 0, axis=1)
    common_phase = (np.max(offsets, axis=1) - np.min(offsets, axis=1)) <= 30.0

    p_order = float(np.sum(weights[ordered]))
    p_common = float(np.sum(weights[common_phase]))
    p_transport = 0.0
    p_coherent_transport = 0.0
    p_reforming_transport = 0.0
    p_ordered_speed_indeterminate = 0.0
    p_incoherent = 0.0
    speed_samples = []
    for jitter in ANCHOR_JITTER_MIN:
        speeds = []
        for i in range(n_nodes - 1):
            dt = activation[:, i + 1] - activation[:, i]
            speeds.append((nodes[i + 1] - nodes[i]) * RSUN_KM / (dt * 60.0))
        dt = float(anchor + jitter) - activation[:, -1]
        speeds.append((3.0 - nodes[-1]) * RSUN_KM / (dt * 60.0))
        speeds = np.stack(speeds, axis=1)
        valid_speed = np.all(
            np.isfinite(speeds)
            & (speeds >= SPEED_RANGE[0])
            & (speeds <= SPEED_RANGE[1]),
            axis=1,
        )
        compatible = ordered & valid_speed
        p_transport += float(np.sum(weights[compatible])) / len(ANCHOR_JITTER_MIN)
        p_coherent_transport += float(
            np.sum(weights[compatible & common_phase])
        ) / len(ANCHOR_JITTER_MIN)
        p_reforming_transport += float(
            np.sum(weights[compatible & ~common_phase])
        ) / len(ANCHOR_JITTER_MIN)
        p_ordered_speed_indeterminate += float(
            np.sum(weights[ordered & ~valid_speed])
        ) / len(ANCHOR_JITTER_MIN)
        p_incoherent += float(np.sum(weights[~ordered])) / len(ANCHOR_JITTER_MIN)
        speed_samples.append(speeds)

    return {
        "p_outward_order": p_order,
        "p_transport_25_300": p_transport,
        "p_common_residual_phase": p_common,
        "scenario_probabilities": {
            "coherent_expanding_transport": p_coherent_transport,
            "reforming_expanding_transport": p_reforming_transport,
            "outward_ordered_speed_indeterminate": p_ordered_speed_indeterminate,
            "timing_incoherent": p_incoherent,
        },
    }


def classify(probabilities):
    p_transport = probabilities["p_transport_25_300"]
    p_order = probabilities["p_outward_order"]
    p_common = probabilities["p_common_residual_phase"]
    if p_transport >= 0.50 and p_common >= 0.50:
        return "coherent_expanding_transport_compatible"
    if p_transport >= 0.50:
        return "moving_or_reforming_expanding_transport_compatible"
    if p_order >= 0.50:
        return "outward_ordered_but_speed_indeterminate"
    return "timing_incoherent"


def analyze_event(offset_products, grid, anchor, model, radii, nodes, lines, area):
    curves = activation_curves(
        offset_products, grid, anchor, model, radii, lines
    )
    product_prob = {
        "total_b": np.zeros_like(curves["total_b"]),
        "pb": np.zeros_like(curves["pb"]),
    }
    joint = np.zeros_like(curves["total_b"])
    overlap = np.zeros(len(nodes), float)
    node_summary = []
    for i, node in enumerate(nodes):
        for product in ("total_b", "pb"):
            base = robust_probability(curves[product][i])
            product_prob[product][i] = broaden_probability(
                base, area["timing_sigma_min"][i]
            )
        overlap[i] = float(np.sum(np.sqrt(
            product_prob["total_b"][i] * product_prob["pb"][i]
        )))
        joint[i] = np.sqrt(
            product_prob["total_b"][i] * product_prob["pb"][i]
        )
        joint[i] /= np.sum(joint[i])
        node_summary.append({
            "radius_rsun": float(node),
            "area_ratio": float(area["area_ratio"][i]),
            "area_source": "measured_pB" if area["measured"][i] else "spherical_fallback",
            "broadening_sigma_min": float(area["timing_sigma_min"][i]),
            "tb_pb_bhattacharyya_overlap": float(overlap[i]),
            **posterior_quantiles(joint[i]),
        })

    tau_nodes = np.interp(nodes, radii, model.tau)
    probabilities = enumerate_sequence_probability(joint, anchor, tau_nodes, nodes)
    result = {
        **probabilities,
        "classification": classify(probabilities),
        "median_tb_pb_overlap": float(np.median(overlap)),
        "nodes": node_summary,
        "area": area,
        "curves": curves,
        "joint_posterior": joint,
        "product_posterior": product_prob,
    }
    return result


def serializable(result):
    return {
        key: value for key, value in result.items()
        if key not in ("area", "curves", "joint_posterior", "product_posterior")
    } | {
        "area": {
            "area_ratio": result["area"]["area_ratio"].tolist(),
            "timing_sigma_min": result["area"]["timing_sigma_min"].tolist(),
            "measured": result["area"]["measured"].tolist(),
            "node_log_area": result["area"]["node_log_area"].tolist(),
        }
    }


def posterior_figure(observed):
    fig, axes = plt.subplots(3, 2, figsize=(13.5, 13.5), constrained_layout=True)
    colors = ("#4c78a8", "#f58518", "#54a24b", "#b279a2")
    for row, representation in enumerate(REPRESENTATIONS):
        for col, candidate_name in enumerate(("event9", "event12")):
            candidate = CANDIDATES[candidate_name]
            result = observed[candidate_name][representation]
            ax = axes[row, col]
            for i, node in enumerate(candidate["nodes"]):
                ax.plot(
                    OFFSETS_MIN, result["joint_posterior"][i], marker="o",
                    color=colors[i], label=f"{node:.3f} R_sun",
                )
            ax.axvline(0, color="black", ls="--", lw=1)
            ax.set(
                title=(f"{representation} | event #{candidate['event_number']} | "
                       f"Ptransport={result['p_transport_25_300']:.2f}"),
                xlabel="Residual activation offset (min)",
                ylabel="Joint tB/pB timing probability",
            )
            ax.legend(fontsize=8, ncol=2)
            ax.grid(alpha=0.2)
    path = OUT / "pds_20080113_expanding_front_timing_posteriors.png"
    fig.savefig(path, dpi=230)
    plt.close(fig)
    return path


def summary_figure(observed, pair_null, primary_p):
    fig, axes = plt.subplots(1, 3, figsize=(16, 5.2), constrained_layout=True)
    labels = [f"#9\n{r}" for r in REPRESENTATIONS] + [f"#12\n{r}" for r in REPRESENTATIONS]
    scenario_keys = (
        "coherent_expanding_transport", "reforming_expanding_transport",
        "outward_ordered_speed_indeterminate", "timing_incoherent",
    )
    scenario_labels = ("coherent transport", "reforming transport",
                       "ordered; speed uncertain", "incoherent")
    scenario_colors = ("#54a24b", "#eeca3b", "#4c78a8", "#e45756")
    bottom = np.zeros(6)
    for key, label, color in zip(scenario_keys, scenario_labels, scenario_colors):
        values = []
        for candidate_name in ("event9", "event12"):
            values.extend([
                observed[candidate_name][r]["scenario_probabilities"][key]
                for r in REPRESENTATIONS
            ])
        axes[0].bar(np.arange(6), values, bottom=bottom, color=color, label=label)
        bottom += np.asarray(values)
    axes[0].set(xticks=np.arange(6), xticklabels=labels, ylim=(0, 1),
                ylabel="Posterior scenario probability",
                title="Expansion-aware timing scenarios")
    axes[0].legend(fontsize=7, loc="lower left")
    axes[0].grid(axis="y", alpha=0.2)

    for candidate_name, marker, color in (("event9", "o", "#4c78a8"),
                                           ("event12", "s", "#f58518")):
        for representation in REPRESENTATIONS:
            nodes = observed[candidate_name][representation]["nodes"]
            axes[1].plot([n["area_ratio"] for n in nodes],
                         [n["radius_rsun"] for n in nodes], marker=marker,
                         color=color, alpha=0.75,
                         label=(f"#{CANDIDATES[candidate_name]['event_number']} {representation}"
                                if representation == "bff" else None))
    axes[1].set(xlabel="Projected area ratio A(r)/A(inner)",
                ylabel="Radius (R_sun)", title="Measured/fallback expansion")
    axes[1].legend(fontsize=8)
    axes[1].grid(alpha=0.2)

    pvals = [pair_null[r]["p_raw"] for r in REPRESENTATIONS]
    axes[2].bar(REPRESENTATIONS, -np.log10(pvals), color="#72b7b2")
    axes[2].axhline(-np.log10(0.05), color="black", ls="--", label="raw p=0.05")
    axes[2].set(ylabel="-log10(raw pair p)",
                title=f"Shifted-pair reference; primary p={primary_p:.3f}")
    axes[2].legend(fontsize=8)
    axes[2].grid(axis="y", alpha=0.2)
    path = OUT / "pds_20080113_expanding_front_summary.png"
    fig.savefig(path, dpi=230)
    plt.close(fig)
    return path


def write_report(summary):
    lines = [
        "# Physics-aware expanding-front sensitivity: events #9 and #12",
        "",
        "## Verdict",
        "",
        summary["verdict"],
        "",
        "This is a post-strict-test, physics-motivated sensitivity analysis. It does not "
        "replace the frozen strict result.",
        "",
        "## Event probabilities (raw shifted-null p; no BH)",
        "",
        "| Event | Representation | P(order) | P(25-300 km/s transport) | "
        "P(common residual phase) | P(coherent/reforming/ordered/incoherent) | "
        "tB/pB overlap | Classification | Individual raw p |",
        "| ---: | --- | ---: | ---: | ---: | ---: | ---: | --- | ---: |",
    ]
    for candidate_name in ("event9", "event12"):
        number = CANDIDATES[candidate_name]["event_number"]
        for representation in REPRESENTATIONS:
            row = summary["events"][candidate_name][representation]
            scenarios = row["scenario_probabilities"]
            lines.append(
                f"| #{number} | {representation} | {row['p_outward_order']:.3f} | "
                f"{row['p_transport_25_300']:.3f} | "
                f"{row['p_common_residual_phase']:.3f} | "
                f"{scenarios['coherent_expanding_transport']:.2f}/"
                f"{scenarios['reforming_expanding_transport']:.2f}/"
                f"{scenarios['outward_ordered_speed_indeterminate']:.2f}/"
                f"{scenarios['timing_incoherent']:.2f} | "
                f"{row['median_tb_pb_overlap']:.3f} | {row['classification']} | "
                f"{row['individual_p_raw']:.4f} |"
            )
    lines += [
        "",
        "## Fixed 360-min shifted-pair reference",
        "",
        "| Representation | Exceedances | Controls | Raw p |",
        "| --- | ---: | ---: | ---: |",
    ]
    for representation in REPRESENTATIONS:
        row = summary["pair_null"][representation]
        lines.append(
            f"| {representation} | {row['exceedances']} | {row['controls']} | "
            f"{row['p_raw']:.4f} |"
        )
    lines += [
        "",
        f"Primary >=2-of-3 exceedances: {summary['primary_pair_exceedances']}/"
        f"{summary['control_pair_count']}; raw p={summary['primary_pair_p_raw']:.4f}.",
        "",
        f"Representations in which both events have P(transport)>=0.50: "
        f"{summary['both_event_compatible_representation_count']}/3.",
        "",
        "## Area-kernel sensitivity at the real events",
        "",
        "| Event | Representation | Measured area Ptransport | Spherical r^2 Ptransport | "
        "Fixed 15-min Ptransport |",
        "| ---: | --- | ---: | ---: | ---: |",
    ]
    for candidate_name in ("event9", "event12"):
        number = CANDIDATES[candidate_name]["event_number"]
        for representation in REPRESENTATIONS:
            measured = summary["events"][candidate_name][representation][
                "p_transport_25_300"
            ]
            spherical = summary["area_kernel_sensitivity"]["spherical"][candidate_name][
                representation
            ]["p_transport_25_300"]
            fixed = summary["area_kernel_sensitivity"]["fixed15"][candidate_name][
                representation
            ]["p_transport_25_300"]
            lines.append(
                f"| #{number} | {representation} | {measured:.3f} | "
                f"{spherical:.3f} | {fixed:.3f} |"
            )
    lines += [
        "",
        f"Maximum absolute change in Ptransport across these area/timing kernels: "
        f"{summary['area_kernel_max_absolute_transport_change']:.3f}.",
        "",
        "## Expansion and node timing",
        "",
    ]
    for candidate_name in ("event9", "event12"):
        number = CANDIDATES[candidate_name]["event_number"]
        lines += [f"### Event #{number}", ""]
        for representation in REPRESENTATIONS:
            lines += [
                f"#### {representation}", "",
                "| R_sun | A/A_inner | Area source | sigma_t (min) | Timing mode | "
                "Mean +/- posterior sigma (min) | 16-84% (min) | tB/pB overlap |",
                "| ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: |",
            ]
            for node in summary["events"][candidate_name][representation]["nodes"]:
                lines.append(
                    f"| {node['radius_rsun']:.3f} | {node['area_ratio']:.3f} | "
                    f"{node['area_source']} | {node['broadening_sigma_min']:.1f} | "
                    f"{node['mode_min']:.0f} | {node['mean_min']:.1f} +/- "
                    f"{node['sigma_min']:.1f} | {node['q16_min']:.0f} to "
                    f"{node['q84_min']:.0f} | "
                    f"{node['tb_pb_bhattacharyya_overlap']:.3f} |"
                )
            lines.append("")
    lines += [
        "## Physical interpretation",
        "",
        summary["physical_interpretation"],
        "",
        "Because every node posterior is normalized separately, radial/spherical brightness "
        "dilution is not penalized. The measured projected area enters only through the "
        "15--30 min timing-broadening kernel. The result therefore tests timing/order "
        "compatibility, not conservation-law closure.",
        "",
        "The analysis does not establish compression, Rankine--Hugoniot consistency, a Mach "
        "number, or an MHD shock branch. Inputs remain Level-0-derived COR1 diagnostics pending "
        "formal SolarSoft SECCHI_PREP Level-1 confirmation.",
        "",
        "All empirical p values are raw; no BH correction is used.",
    ]
    path = OUT / "PDS_20080113_expanding_front_sensitivity_no_BH.md"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def main():
    grid, radii, event_times, products = strict.prepare_inputs()
    models = {
        name: phase.model_from_dict(candidate["model"], radii)
        for name, candidate in CANDIDATES.items()
    }
    lines, _ = strict.observed_geometries(products, grid, radii, event_times, models)

    with np.load(geometry.JAN_MAPS) as data:
        path_pa = np.asarray(data["traced_path_pa"], float)
    moment = geometry.streamer_moments(geometry.JAN_CUBE, "pb_signed", path_pa)

    observed = {name: {} for name in CANDIDATES}
    observed_area = {}
    for candidate_name, candidate in CANDIDATES.items():
        anchor = float(event_times[candidate["event_index"]])
        area = event_area(moment, anchor, models[candidate_name], radii, candidate["nodes"])
        observed_area[candidate_name] = area
        for representation in REPRESENTATIONS:
            observed[candidate_name][representation] = analyze_event(
                products[representation], grid, anchor, models[candidate_name], radii,
                candidate["nodes"], lines[candidate_name][representation], area,
            )

    area_kernel_sensitivity = {"spherical": {}, "fixed15": {}}
    max_area_kernel_change = 0.0
    for scheme in area_kernel_sensitivity:
        for candidate_name, candidate in CANDIDATES.items():
            nodes = candidate["nodes"]
            anchor = float(event_times[candidate["event_index"]])
            if scheme == "spherical":
                ratio = np.square(nodes / nodes[0])
                sigma = np.clip(15.0 * np.sqrt(ratio), 15.0, 30.0)
            else:
                ratio = np.ones(len(nodes), float)
                sigma = np.full(len(nodes), 15.0)
            area = {
                "area_ratio": ratio,
                "timing_sigma_min": sigma,
                "measured": np.zeros(len(nodes), bool),
                "node_log_area": np.full(len(nodes), np.nan),
            }
            area_kernel_sensitivity[scheme][candidate_name] = {}
            for representation in REPRESENTATIONS:
                result = analyze_event(
                    products[representation], grid, anchor, models[candidate_name], radii,
                    nodes, lines[candidate_name][representation], area,
                )
                values = {
                    "p_outward_order": result["p_outward_order"],
                    "p_transport_25_300": result["p_transport_25_300"],
                    "p_common_residual_phase": result["p_common_residual_phase"],
                }
                area_kernel_sensitivity[scheme][candidate_name][representation] = values
                max_area_kernel_change = max(
                    max_area_kernel_change,
                    abs(values["p_transport_25_300"]
                        - observed[candidate_name][representation]["p_transport_25_300"]),
                )

    with PAIR_NULL_FILE.open(newline="", encoding="utf-8") as handle:
        early_anchors = np.array([
            float(row["early_anchor_min"]) for row in csv.DictReader(handle)
        ])

    null_rows = []
    for number, early_anchor in enumerate(early_anchors, 1):
        late_anchor = float(early_anchor + strict.PAIR_SEPARATION_MIN)
        areas = {
            "event9": event_area(
                moment, early_anchor, models["event9"], radii,
                CANDIDATES["event9"]["nodes"],
            ),
            "event12": event_area(
                moment, late_anchor, models["event12"], radii,
                CANDIDATES["event12"]["nodes"],
            ),
        }
        row = {
            "early_anchor_min": early_anchor,
            "late_anchor_min": late_anchor,
            "early_utc": utc_from_minute(early_anchor),
            "late_utc": utc_from_minute(late_anchor),
        }
        count = 0
        for representation in REPRESENTATIONS:
            early = analyze_event(
                products[representation], grid, early_anchor, models["event9"], radii,
                CANDIDATES["event9"]["nodes"], lines["event9"][representation],
                areas["event9"],
            )
            late = analyze_event(
                products[representation], grid, late_anchor, models["event12"], radii,
                CANDIDATES["event12"]["nodes"], lines["event12"][representation],
                areas["event12"],
            )
            early_score = early["p_transport_25_300"]
            late_score = late["p_transport_25_300"]
            passed = bool(
                early_score >= observed["event9"][representation]["p_transport_25_300"]
                and late_score >= observed["event12"][representation]["p_transport_25_300"]
            )
            row[f"{representation}_early_p_transport"] = early_score
            row[f"{representation}_late_p_transport"] = late_score
            row[f"{representation}_pair_exceed"] = passed
            count += int(passed)
        row["representation_exceedance_count"] = count
        row["primary_pair_exceed"] = bool(count >= 2)
        null_rows.append(row)
        if number % 20 == 0:
            print(f"Expanding-front pair null: {number}/{len(early_anchors)}", flush=True)

    pair_null = {}
    event_p = {name: {} for name in CANDIDATES}
    for representation in REPRESENTATIONS:
        pair_k = int(sum(row[f"{representation}_pair_exceed"] for row in null_rows))
        early_k = int(sum(
            row[f"{representation}_early_p_transport"]
            >= observed["event9"][representation]["p_transport_25_300"]
            for row in null_rows
        ))
        late_k = int(sum(
            row[f"{representation}_late_p_transport"]
            >= observed["event12"][representation]["p_transport_25_300"]
            for row in null_rows
        ))
        pair_null[representation] = {
            "exceedances": pair_k,
            "controls": len(null_rows),
            "p_raw": empirical_p(pair_k, len(null_rows)),
        }
        event_p["event9"][representation] = empirical_p(early_k, len(null_rows))
        event_p["event12"][representation] = empirical_p(late_k, len(null_rows))

    primary_k = int(sum(row["primary_pair_exceed"] for row in null_rows))
    primary_p = empirical_p(primary_k, len(null_rows))
    compatible_reps = int(sum(
        observed["event9"][rep]["p_transport_25_300"] >= 0.50
        and observed["event12"][rep]["p_transport_25_300"] >= 0.50
        for rep in REPRESENTATIONS
    ))
    passed = bool(primary_p < 0.05 and compatible_reps >= 2)

    if passed:
        verdict = (
            "Under the declared expansion-aware timing model, the ordered #9 -> #12 pair "
            "is processing-robustly compatible with outward 25--300 km/s pattern transport "
            "in at least two representations and is locally unusual in the shifted-pair "
            "reference. This rescues transport compatibility, not a shock identification."
        )
    elif primary_p < 0.05:
        verdict = (
            "The expansion-aware timing score of the ordered #9 -> #12 pair is locally "
            "unusual in the >=2-of-3 shifted reference (raw p<0.05), with independent support "
            "from base60 and nrgf60. However, neither event pair reaches the predeclared "
            "P(25--300 km/s transport)>=0.50 threshold in two representations. The positive "
            "result is therefore an outward-ordered, phase-jittered morphology indication, "
            "not a completed transport or shock detection."
        )
    else:
        verdict = (
            "Accounting for projected expansion and timing broadening does not satisfy the "
            "complete pair rule. The morphology may remain compatible with individual "
            "outward sequences, but the pair is not a processing-robust transport detection."
        )

    classifications = {
        name: [observed[name][rep]["classification"] for rep in REPRESENTATIONS]
        for name in CANDIDATES
    }
    reforming = any("reforming" in value for values in classifications.values() for value in values)
    if passed and reforming:
        physical_interpretation = (
            "The posterior favors an outward but non-rigid sequence: node timing is compatible "
            "with transport after expansion broadening, while residual phases are not uniformly "
            "locked. The appropriate phenomenology is therefore a moving/reforming cusp-gate "
            "or shock-cell-like pattern, not one compact ballistic pulse and not a universal "
            "stationary Habbal shock."
        )
    elif passed:
        physical_interpretation = (
            "The posterior is compatible with a common expanding propagation sequence across "
            "the frozen nodes. It remains an imaging-only pattern result rather than a physical "
            "shock closure."
        )
    elif primary_p < 0.05:
        physical_interpretation = (
            "The two filter-robust representations favor an outward-ordered but phase-jittered "
            "sequence. Event #9 is dominated by correct radial ordering with poorly constrained "
            "inter-node speeds; event #12 lies close to full expanding-transport compatibility. "
            "The combined behavior is more consistent with a moving/reforming cusp-gate than "
            "with one rigid ballistic pulse or a universal stationary Habbal shock. The evidence "
            "is suggestive rather than a closed causal chain."
        )
    else:
        physical_interpretation = (
            "The sensitivity does not close the transport chain. The safe interpretation "
            "remains processing-robust X/diamond-like morphology with event-dependent timing, "
            "consistent with a moving/reforming gate but not uniquely distinguishable from "
            "line-of-sight/topological structure."
        )

    event_summary = {name: {} for name in CANDIDATES}
    for candidate_name in CANDIDATES:
        for representation in REPRESENTATIONS:
            row = serializable(observed[candidate_name][representation])
            row["individual_p_raw"] = event_p[candidate_name][representation]
            event_summary[candidate_name][representation] = row

    summary = {
        "analysis_status": "post-strict-test physics-motivated sensitivity",
        "events": event_summary,
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
        "offsets_min": OFFSETS_MIN.tolist(),
        "anchor_jitter_min": ANCHOR_JITTER_MIN.tolist(),
        "speed_range_km_s": list(SPEED_RANGE),
        "pair_null": pair_null,
        "area_kernel_sensitivity": area_kernel_sensitivity,
        "area_kernel_max_absolute_transport_change": max_area_kernel_change,
        "control_pair_count": len(null_rows),
        "primary_pair_exceedances": primary_k,
        "primary_pair_p_raw": primary_p,
        "both_event_compatible_representation_count": compatible_reps,
        "passed_expanding_transport_compatibility_rule": passed,
        "verdict": verdict,
        "physical_interpretation": physical_interpretation,
        "no_BH": True,
        "calibration": "Level-0-derived X response plus calibrated-port pB width diagnostic",
    }

    post_figure = posterior_figure(observed)
    summary_plot = summary_figure(observed, pair_null, primary_p)
    report_path = write_report(summary)
    summary["files"] = {
        "posterior_figure": post_figure.name,
        "summary_figure": summary_plot.name,
        "report": report_path.name,
    }
    json_path = OUT / "pds_20080113_expanding_front_sensitivity_results_no_BH.json"
    json_path.write_text(json.dumps(summary, indent=2, allow_nan=True), encoding="utf-8")

    csv_path = OUT / "pds_20080113_expanding_front_pair_null.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(null_rows[0]))
        writer.writeheader()
        writer.writerows(null_rows)

    payload = {"offsets_min": OFFSETS_MIN}
    for candidate_name in CANDIDATES:
        for representation in REPRESENTATIONS:
            result = observed[candidate_name][representation]
            payload[f"{candidate_name}_{representation}_joint_posterior"] = result[
                "joint_posterior"
            ]
            payload[f"{candidate_name}_{representation}_tb_curve"] = result["curves"][
                "total_b"
            ]
            payload[f"{candidate_name}_{representation}_pb_curve"] = result["curves"]["pb"]
    np.savez_compressed(
        OUT / "pds_20080113_expanding_front_posteriors.npz", **payload
    )

    print(json.dumps({
        "passed": passed,
        "primary_pair_p_raw": primary_p,
        "compatible_representations": compatible_reps,
        "events": {
            name: {
                rep: {
                    "p_order": observed[name][rep]["p_outward_order"],
                    "p_transport": observed[name][rep]["p_transport_25_300"],
                    "p_common": observed[name][rep]["p_common_residual_phase"],
                    "classification": observed[name][rep]["classification"],
                    "individual_p_raw": event_p[name][rep],
                }
                for rep in REPRESENTATIONS
            }
            for name in CANDIDATES
        },
        "pair_null": pair_null,
    }, indent=2))


if __name__ == "__main__":
    main()
