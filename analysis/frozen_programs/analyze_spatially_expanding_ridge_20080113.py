#!/usr/bin/env python3
"""Spatially expanding-ridge pilot for COR1-A events #9 and #12.

The frozen method is documented in
PDS_20080113_spatially_expanding_ridge_pilot_spec.md.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import analyze_dynamic_geometry_gate as geometry
import analyze_expanding_front_sensitivity_20080113 as timing
import analyze_event09_2d_xfront_20080113 as x2d
import analyze_node_activation_order_20080113 as strict
import analyze_pds_event_phase_jitter_20080111_14 as phase


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "results" / "pds_20080111_14_nonlinear_level0"
OUT = SOURCE / "spatially_expanding_ridge_pilot"
OUT.mkdir(parents=True, exist_ok=True)

SCHEMES = ("fixed", "spherical", "measured")
REPRESENTATIONS = strict.REPRESENTATIONS
CANDIDATES = strict.CANDIDATES
OFFSETS_MIN = timing.OFFSETS_MIN
PAIR_NULL_FILE = strict.PAIR_NULL_FILE
BASE_SIGMA_RSUN = 0.025
MAX_SIGMA_RSUN = 0.075
GAUSSIAN_U = np.linspace(-3.0, 3.0, 13)
GAUSSIAN_W = np.exp(-0.5 * GAUSSIAN_U**2)


def json_default(value):
    if isinstance(value, np.bool_):
        return bool(value)
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.floating):
        return float(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    raise TypeError(f"Not JSON serializable: {type(value).__name__}")


def scheme_area(moment, anchor, model, radii, nodes, scheme):
    nodes = np.asarray(nodes, float)
    if scheme == "fixed":
        ratio = np.ones(len(nodes), float)
        sigma_t = np.full(len(nodes), 15.0)
        measured = np.zeros(len(nodes), bool)
    elif scheme == "spherical":
        ratio = np.square(nodes / nodes[0])
        sigma_t = np.clip(15.0 * np.sqrt(ratio), 15.0, 30.0)
        measured = np.zeros(len(nodes), bool)
    elif scheme == "measured":
        row = timing.event_area(moment, anchor, model, radii, nodes)
        ratio = np.asarray(row["area_ratio"], float)
        sigma_t = np.asarray(row["timing_sigma_min"], float)
        measured = np.asarray(row["measured"], bool)
    else:
        raise ValueError(scheme)
    return {
        "scheme": scheme,
        "nodes": nodes,
        "area_ratio": ratio,
        "timing_sigma_min": sigma_t,
        "measured": measured,
    }


def area_ratio_at(radius, area):
    radius = np.asarray(radius, float)
    scheme = area["scheme"]
    if scheme == "fixed":
        return np.ones_like(radius)
    if scheme == "spherical":
        return np.square(radius / area["nodes"][0])
    ratio = np.interp(
        radius,
        area["nodes"],
        area["area_ratio"],
        left=area["area_ratio"][0],
        right=area["area_ratio"][-1],
    )
    spherical = np.square(radius / area["nodes"][0])
    return np.where(np.isfinite(ratio) & (ratio > 0), ratio, spherical)


def sigma_r_at(radius, area):
    ratio = np.maximum(area_ratio_at(radius, area), 1e-6)
    return np.clip(BASE_SIGMA_RSUN * np.sqrt(ratio), BASE_SIGMA_RSUN, MAX_SIGMA_RSUN)


def radial_sample(map2d, radii, targets, pa_indices):
    """Vectorized linear interpolation along radius at fixed PA indices."""
    targets = np.asarray(targets, float)
    pidx = np.broadcast_to(np.asarray(pa_indices, int), targets.shape)
    step = float(np.nanmedian(np.diff(radii)))
    position = (targets - radii[0]) / step
    valid = np.isfinite(position) & (position >= 0) & (position <= len(radii) - 1)
    lo_raw = np.floor(np.where(valid, position, 0.0)).astype(int)
    lo = np.clip(lo_raw, 0, len(radii) - 1)
    hi = np.clip(lo + 1, 0, len(radii) - 1)
    weight = position - lo_raw
    values = (1.0 - weight) * map2d[lo, pidx] + weight * map2d[hi, pidx]
    values[~valid] = np.nan
    return values


def gaussian_radial_mean(map2d, radii, centers, sigma, pa_indices):
    centers = np.asarray(centers, float)
    sigma = np.asarray(sigma, float)
    targets = centers[:, None] + sigma[:, None] * GAUSSIAN_U[None, :]
    pidx = np.broadcast_to(np.asarray(pa_indices, int)[:, None], targets.shape)
    values = radial_sample(map2d, radii, targets, pidx)
    weights = np.broadcast_to(GAUSSIAN_W[None, :], values.shape).copy()
    weights[~np.isfinite(values)] = 0.0
    denominator = np.sum(weights, axis=1)
    numerator = np.nansum(values * weights, axis=1)
    output = np.full(len(centers), np.nan)
    good = denominator > 0
    output[good] = numerator[good] / denominator[good]
    return output


def expanding_diagonal_response(map2d, radii, vertex, signed_slope, area):
    pa_indices = np.flatnonzero(x2d.ARM_MASK)
    q = x2d.PA_OFFSETS[pa_indices]
    target = float(vertex) + float(signed_slope) * q
    sigma = sigma_r_at(target, area)
    side_distance = np.maximum(x2d.SIDEBAND_RSUN, 3.0 * sigma)
    line = gaussian_radial_mean(map2d, radii, target, sigma, pa_indices)
    side_plus = gaussian_radial_mean(
        map2d, radii, target + side_distance, sigma, pa_indices
    )
    side_minus = gaussian_radial_mean(
        map2d, radii, target - side_distance, sigma, pa_indices
    )
    valid = np.isfinite(line) & np.isfinite(side_plus) & np.isfinite(side_minus)
    if np.sum(valid) < x2d.MIN_ARM_SAMPLES:
        return np.nan
    contrast = line[valid] - 0.5 * (side_plus[valid] + side_minus[valid])
    return float(np.mean(contrast))


def expanding_x_response(maps, radii, geometry_row, area):
    vertex = float(geometry_row["vertex_radius_rsun"])
    slope = float(geometry_row["absolute_slope_rsun_per_degree"])
    values = {}
    for product, map2d in maps.items():
        positive = expanding_diagonal_response(map2d, radii, vertex, slope, area)
        negative = expanding_diagonal_response(map2d, radii, vertex, -slope, area)
        values[product] = min(positive, negative) if (
            np.isfinite(positive) and np.isfinite(negative)
        ) else np.nan
    coherent = min(values.values()) if all(np.isfinite(list(values.values()))) else np.nan
    return values, coherent


def fit_node_geometry(maps, radii, nominal_node, area):
    best = None
    for vertex in float(nominal_node) + x2d.VERTEX_OFFSETS:
        for slope in x2d.SLOPES:
            geometry_row = {
                "nominal_radius_rsun": float(nominal_node),
                "vertex_radius_rsun": float(vertex),
                "absolute_slope_rsun_per_degree": float(slope),
            }
            values, coherent = expanding_x_response(maps, radii, geometry_row, area)
            candidate = {
                **geometry_row,
                "tb_x_response": float(values["total_b"]),
                "pb_x_response": float(values["pb"]),
                "coherent_x_response": float(coherent),
                "sigma_r_rsun": float(sigma_r_at(np.array([nominal_node]), area)[0]),
            }
            if best is None or (
                np.isfinite(coherent)
                and (
                    not np.isfinite(best["coherent_x_response"])
                    or coherent > best["coherent_x_response"]
                )
            ):
                best = candidate
    return best


def fit_event_geometry(products, grid, anchor, model, radii, nodes, area):
    maps = {
        product: x2d.aligned_slice(cube, grid, float(anchor), model.tau)
        for product, cube in products.items()
    }
    return [fit_node_geometry(maps, radii, node, area) for node in nodes], maps


def activation_curves(products, grid, anchor, model, radii, geometries, area):
    curves = {
        "total_b": np.full((len(geometries), len(OFFSETS_MIN)), np.nan),
        "pb": np.full((len(geometries), len(OFFSETS_MIN)), np.nan),
    }
    for k, offset in enumerate(OFFSETS_MIN):
        maps = {
            product: x2d.aligned_slice(cube, grid, float(anchor + offset), model.tau)
            for product, cube in products.items()
        }
        for i, geometry_row in enumerate(geometries):
            values, _ = expanding_x_response(maps, radii, geometry_row, area)
            curves["total_b"][i, k] = values["total_b"]
            curves["pb"][i, k] = values["pb"]
    return curves


def safe_curve_probability(curve):
    """Return an uninformative posterior when a broadened control is unsampled."""
    probability = timing.robust_probability(curve)
    total = float(np.nansum(probability))
    if not np.isfinite(total) or total <= 0:
        return np.full(len(OFFSETS_MIN), 1.0 / len(OFFSETS_MIN))
    probability = np.where(np.isfinite(probability), probability, 0.0)
    total = float(np.sum(probability))
    if total <= 0:
        return np.full(len(OFFSETS_MIN), 1.0 / len(OFFSETS_MIN))
    return probability / total


def analyze_event(products, grid, anchor, model, radii, nodes, geometries, area):
    curves = activation_curves(
        products, grid, anchor, model, radii, geometries, area
    )
    product_p = {name: np.zeros_like(curves[name]) for name in ("total_b", "pb")}
    joint = np.zeros_like(curves["total_b"])
    overlaps = np.zeros(len(nodes), float)
    node_rows = []
    for i, node in enumerate(nodes):
        for product in ("total_b", "pb"):
            base = safe_curve_probability(curves[product][i])
            product_p[product][i] = timing.broaden_probability(
                base, area["timing_sigma_min"][i]
            )
        overlaps[i] = float(np.sum(np.sqrt(product_p["total_b"][i] * product_p["pb"][i])))
        joint[i] = np.sqrt(product_p["total_b"][i] * product_p["pb"][i])
        joint[i] /= np.sum(joint[i])
        node_rows.append({
            "radius_rsun": float(node),
            "area_ratio": float(area["area_ratio"][i]),
            "timing_sigma_min": float(area["timing_sigma_min"][i]),
            "spatial_sigma_rsun": float(sigma_r_at(np.array([node]), area)[0]),
            "tb_pb_overlap": float(overlaps[i]),
            **timing.posterior_quantiles(joint[i]),
        })
    tau_nodes = np.interp(nodes, radii, model.tau)
    probabilities = timing.enumerate_sequence_probability(
        joint, anchor, tau_nodes, nodes
    )
    zero = int(np.argmin(np.abs(OFFSETS_MIN)))
    median_zero = float(np.nanmedian(np.minimum(
        curves["total_b"][:, zero], curves["pb"][:, zero]
    )))
    return {
        **probabilities,
        "classification": timing.classify(probabilities),
        "median_tb_pb_overlap": float(np.nanmedian(overlaps)),
        "median_zero_offset_x": median_zero,
        "nodes": node_rows,
        "geometries": geometries,
        "curves": curves,
        "joint_posterior": joint,
    }


def serializable_event(row):
    return {key: value for key, value in row.items() if key not in ("curves", "joint_posterior")}


def coherent_map(maps):
    same = maps["total_b"] * maps["pb"] > 0
    output = np.zeros_like(maps["total_b"])
    output[same] = np.sign(maps["total_b"][same]) * np.sqrt(
        np.abs(maps["total_b"][same] * maps["pb"][same])
    )
    return output


def make_map_figure(observed_maps, observed, radii):
    fig, axes = plt.subplots(3, 2, figsize=(13, 14), constrained_layout=True)
    for row_index, scheme in enumerate(SCHEMES):
        for col_index, candidate_name in enumerate(("event9", "event12")):
            ax = axes[row_index, col_index]
            maps = observed_maps[candidate_name]["base60"]
            image = ax.pcolormesh(
                x2d.PA_OFFSETS, radii, coherent_map(maps), cmap="RdBu_r",
                vmin=-2.5, vmax=2.5, shading="auto",
            )
            result = observed[scheme][candidate_name]["base60"]
            area = result["_area"]
            for geometry_row in result["geometries"]:
                q = np.linspace(-8.0, 8.0, 121)
                vertex = geometry_row["vertex_radius_rsun"]
                slope = geometry_row["absolute_slope_rsun_per_degree"]
                for sign, color in ((1.0, "#f6bd60"), (-1.0, "#84a59d")):
                    center = vertex + sign * slope * q
                    width = sigma_r_at(center, area)
                    ax.plot(q, center, color=color, lw=1.1)
                    ax.fill_between(q, center - width, center + width,
                                    color=color, alpha=0.12, linewidth=0)
            ax.set(
                title=(f"{scheme} | #{CANDIDATES[candidate_name]['event_number']} | "
                       f"Ptransport={result['p_transport_25_300']:.3f}"),
                xlabel="PA offset (deg)", ylabel=r"Radius ($R_\odot$)",
            )
            fig.colorbar(image, ax=ax, label="Same-sign tB/pB z")
    path = OUT / "pds_20080113_spatially_expanding_ridge_maps.png"
    fig.savefig(path, dpi=220)
    plt.close(fig)
    return path


def make_summary_figure(observed, pair_summary):
    fig, axes = plt.subplots(1, 3, figsize=(16, 5.5), constrained_layout=True)
    x = np.arange(len(REPRESENTATIONS))
    width = 0.25
    colors = {"fixed": "#4c78a8", "spherical": "#f2cf5b", "measured": "#e45756"}
    for ax, candidate_name in zip(axes[:2], ("event9", "event12")):
        for j, scheme in enumerate(SCHEMES):
            values = [observed[scheme][candidate_name][rep]["p_transport_25_300"]
                      for rep in REPRESENTATIONS]
            ax.bar(x + (j - 1) * width, values, width, label=scheme,
                   color=colors[scheme])
        ax.axhline(0.5, color="black", ls="--", lw=1)
        ax.set(
            xticks=x, xticklabels=REPRESENTATIONS, ylim=(0, 1),
            ylabel="P(25-300 km/s transport)",
            title=f"Event #{CANDIDATES[candidate_name]['event_number']}",
        )
        ax.grid(axis="y", alpha=0.2)
    axes[0].legend(fontsize=8)
    pvalues = [pair_summary[scheme]["primary_pair_p_raw"] for scheme in SCHEMES]
    axes[2].bar(SCHEMES, -np.log10(pvalues), color=[colors[s] for s in SCHEMES])
    axes[2].axhline(-np.log10(0.05), color="black", ls="--", lw=1,
                    label="raw p=0.05")
    axes[2].set(ylabel="-log10(raw pair p)", title="131 shifted pairs")
    axes[2].legend(fontsize=8)
    axes[2].grid(axis="y", alpha=0.2)
    path = OUT / "pds_20080113_spatially_expanding_ridge_summary.png"
    fig.savefig(path, dpi=220)
    plt.close(fig)
    return path


def write_report(summary):
    lines = [
        "# Spatially expanding-ridge pilot: events #9 and #12",
        "",
        "## Verdict",
        "",
        summary["verdict"],
        "",
        "This is a frozen post-extraction sensitivity pilot.  Fixed, spherical, and "
        "measured-width filters are a sensitivity bracket, not independent discoveries.",
        "",
        "## Event probabilities",
        "",
        "| Width law | Event | Representation | P(order) | P(transport) | P(common phase) | "
        "Median spatial sigma (R_sun) | Classification |",
        "| --- | ---: | --- | ---: | ---: | ---: | ---: | --- |",
    ]
    for scheme in SCHEMES:
        for candidate_name in ("event9", "event12"):
            for representation in REPRESENTATIONS:
                row = summary["events"][scheme][candidate_name][representation]
                spatial = np.median([node["spatial_sigma_rsun"] for node in row["nodes"]])
                lines.append(
                    f"| {scheme} | #{CANDIDATES[candidate_name]['event_number']} | "
                    f"{representation} | {row['p_outward_order']:.3f} | "
                    f"{row['p_transport_25_300']:.3f} | "
                    f"{row['p_common_residual_phase']:.3f} | {spatial:.3f} | "
                    f"{row['classification']} |"
                )
    lines += [
        "",
        "## Fixed 360-min pair reference (raw; no BH)",
        "",
        "| Width law | Representation | Exceedances / 131 | Raw p | Primary >=2/3 |",
        "| --- | --- | ---: | ---: | ---: |",
    ]
    for scheme in SCHEMES:
        for representation in REPRESENTATIONS:
            row = summary["pair_null"][scheme]["representations"][representation]
            lines.append(
                f"| {scheme} | {representation} | {row['exceedances']}/131 | "
                f"{row['p_raw']:.4f} |  |"
            )
        row = summary["pair_null"][scheme]
        lines.append(
            f"| {scheme} | primary | {row['primary_pair_exceedances']}/131 | "
            f"{row['primary_pair_p_raw']:.4f} | yes |"
        )
    lines += [
        "",
        "## Geometry changes relative to the one-pixel fixed-width filter",
        "",
        "| Event | Representation | Width law | Maximum |delta vertex| (R_sun) | "
        "Changed slopes / nodes |",
        "| ---: | --- | --- | ---: | ---: |",
    ]
    for candidate_name in ("event9", "event12"):
        for representation in REPRESENTATIONS:
            reference = summary["events"]["fixed"][candidate_name][representation]["geometries"]
            for scheme in ("spherical", "measured"):
                current = summary["events"][scheme][candidate_name][representation]["geometries"]
                delta = max(abs(a["vertex_radius_rsun"] - b["vertex_radius_rsun"])
                            for a, b in zip(reference, current))
                changed = sum(a["absolute_slope_rsun_per_degree"] !=
                              b["absolute_slope_rsun_per_degree"]
                              for a, b in zip(reference, current))
                lines.append(
                    f"| #{CANDIDATES[candidate_name]['event_number']} | {representation} | "
                    f"{scheme} | {delta:.3f} | {changed}/{len(current)} |"
                )
    lines += [
        "",
        "## Interpretation limits",
        "",
        summary["physical_interpretation"],
        "",
        "The pB projected width [r sigma_PA]^2 is an imaging proxy rather than a magnetic "
        "flux-tube area measurement.  The Gaussian radial-normal approximation is therefore "
        "a controlled sensitivity model, not a unique MHD forward model.",
        "",
        "No density jump, Rankine--Hugoniot closure, magnetosonic Mach number, or MHD shock "
        "branch is inferred.  All empirical p values are raw; no BH correction is used.",
    ]
    path = OUT / "PDS_20080113_spatially_expanding_ridge_pilot_no_BH.md"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def main():
    grid, radii, event_times, products = strict.prepare_inputs()
    models = {
        name: phase.model_from_dict(candidate["model"], radii)
        for name, candidate in CANDIDATES.items()
    }
    with np.load(geometry.JAN_MAPS) as data:
        path_pa = np.asarray(data["traced_path_pa"], float)
    moment = geometry.streamer_moments(geometry.JAN_CUBE, "pb_signed", path_pa)

    observed = {scheme: {name: {} for name in CANDIDATES} for scheme in SCHEMES}
    observed_maps = {name: {} for name in CANDIDATES}
    fitted_geometries = {scheme: {name: {} for name in CANDIDATES} for scheme in SCHEMES}
    observed_areas = {scheme: {} for scheme in SCHEMES}

    for scheme in SCHEMES:
        for candidate_name, candidate in CANDIDATES.items():
            anchor = float(event_times[candidate["event_index"]])
            area = scheme_area(
                moment, anchor, models[candidate_name], radii, candidate["nodes"], scheme
            )
            observed_areas[scheme][candidate_name] = area
            for representation in REPRESENTATIONS:
                geometries, maps = fit_event_geometry(
                    products[representation], grid, anchor, models[candidate_name],
                    radii, candidate["nodes"], area,
                )
                fitted_geometries[scheme][candidate_name][representation] = geometries
                if scheme == "fixed":
                    observed_maps[candidate_name][representation] = maps
                result = analyze_event(
                    products[representation], grid, anchor, models[candidate_name],
                    radii, candidate["nodes"], geometries, area,
                )
                result["_area"] = area
                observed[scheme][candidate_name][representation] = result

    with PAIR_NULL_FILE.open(newline="", encoding="utf-8") as handle:
        early_anchors = np.array([
            float(row["early_anchor_min"]) for row in csv.DictReader(handle)
        ])

    null_rows = []
    for index, early_anchor in enumerate(early_anchors, 1):
        late_anchor = float(early_anchor + strict.PAIR_SEPARATION_MIN)
        row = {"early_anchor_min": float(early_anchor), "late_anchor_min": late_anchor}
        for scheme in SCHEMES:
            control_areas = {
                "event9": scheme_area(
                    moment, early_anchor, models["event9"], radii,
                    CANDIDATES["event9"]["nodes"], scheme,
                ),
                "event12": scheme_area(
                    moment, late_anchor, models["event12"], radii,
                    CANDIDATES["event12"]["nodes"], scheme,
                ),
            }
            count = 0
            for representation in REPRESENTATIONS:
                early = analyze_event(
                    products[representation], grid, early_anchor, models["event9"], radii,
                    CANDIDATES["event9"]["nodes"],
                    fitted_geometries[scheme]["event9"][representation],
                    control_areas["event9"],
                )
                late = analyze_event(
                    products[representation], grid, late_anchor, models["event12"], radii,
                    CANDIDATES["event12"]["nodes"],
                    fitted_geometries[scheme]["event12"][representation],
                    control_areas["event12"],
                )
                early_score = early["p_transport_25_300"]
                late_score = late["p_transport_25_300"]
                passed = bool(
                    early_score >= observed[scheme]["event9"][representation]["p_transport_25_300"]
                    and late_score >= observed[scheme]["event12"][representation]["p_transport_25_300"]
                )
                prefix = f"{scheme}_{representation}"
                row[f"{prefix}_early_p_transport"] = early_score
                row[f"{prefix}_late_p_transport"] = late_score
                row[f"{prefix}_pair_exceed"] = passed
                count += int(passed)
            row[f"{scheme}_representation_exceedance_count"] = count
            row[f"{scheme}_primary_pair_exceed"] = bool(count >= 2)
        null_rows.append(row)
        if index % 20 == 0:
            print(f"Spatial ridge pair null: {index}/{len(early_anchors)}", flush=True)

    pair_summary = {}
    for scheme in SCHEMES:
        representations = {}
        for representation in REPRESENTATIONS:
            key = f"{scheme}_{representation}_pair_exceed"
            count = int(sum(bool(row[key]) for row in null_rows))
            representations[representation] = {
                "exceedances": count,
                "controls": len(null_rows),
                "p_raw": timing.empirical_p(count, len(null_rows)),
            }
        primary_count = int(sum(
            bool(row[f"{scheme}_primary_pair_exceed"]) for row in null_rows
        ))
        pair_summary[scheme] = {
            "representations": representations,
            "primary_pair_exceedances": primary_count,
            "primary_pair_p_raw": timing.empirical_p(primary_count, len(null_rows)),
        }

    measured_stable = sum(
        observed["measured"]["event9"][rep]["p_outward_order"] >= 0.50
        and observed["measured"]["event12"][rep]["p_outward_order"] >= 0.50
        for rep in REPRESENTATIONS
    ) >= 2
    measured_transport = sum(
        observed["measured"]["event9"][rep]["p_transport_25_300"] >= 0.50
        and observed["measured"]["event12"][rep]["p_transport_25_300"] >= 0.50
        for rep in REPRESENTATIONS
    ) >= 2
    measured_unusual = pair_summary["measured"]["primary_pair_p_raw"] < 0.05

    if measured_transport and measured_unusual:
        verdict = (
            "The measured-width spatial filter makes the #9 -> #12 pair compatible with "
            "25--300 km/s outward transport in at least two representations and unusual "
            "against the frozen 131-pair reference.  Spatial expansion therefore strengthens "
            "transport compatibility, but does not identify an MHD shock."
        )
    elif measured_stable and measured_unusual:
        verdict = (
            "The measured-width spatial filter preserves an unusual outward-ordered #9 -> #12 "
            "pair, but the full transport probability does not cross 0.50 for both events in "
            "two representations.  Spatial broadening supports a moving/reforming morphology, "
            "not a closed ballistic or shock chain."
        )
    elif measured_stable:
        verdict = (
            "Spatial expansion preserves outward ordering, but the measured-width pair is not "
            "unusual in the frozen shifted reference.  The morphology remains compatible but "
            "not statistically distinctive."
        )
    else:
        verdict = (
            "Allowing the ridge to broaden spatially does not preserve a two-representation "
            "ordered pair.  The safe result remains an event-dependent morphology indication."
        )

    max_geometry_shift = 0.0
    for candidate_name in CANDIDATES:
        for representation in REPRESENTATIONS:
            fixed_rows = fitted_geometries["fixed"][candidate_name][representation]
            measured_rows = fitted_geometries["measured"][candidate_name][representation]
            max_geometry_shift = max(
                max_geometry_shift,
                max(abs(a["vertex_radius_rsun"] - b["vertex_radius_rsun"])
                    for a, b in zip(fixed_rows, measured_rows)),
            )

    physical_interpretation = (
        "A result stable from the one-pixel filter through spherical and measured broadening "
        "is not an artifact of imposing a rigid ridge thickness.  A change confined to the "
        "measured filter instead indicates sensitivity to streamer geometry.  In either case, "
        "phase-jittered ordered morphology is most naturally described as an expanding or "
        "reforming cusp/current-sheet gate until density-jump and Mach-number measurements exist."
    )

    public_events = {scheme: {name: {} for name in CANDIDATES} for scheme in SCHEMES}
    for scheme in SCHEMES:
        for candidate_name in CANDIDATES:
            for representation in REPRESENTATIONS:
                public_events[scheme][candidate_name][representation] = serializable_event(
                    observed[scheme][candidate_name][representation]
                )

    previous_path = SOURCE / "expanding_front_sensitivity" / \
        "pds_20080113_expanding_front_sensitivity_results_no_BH.json"
    previous = json.loads(previous_path.read_text(encoding="utf-8"))
    summary = {
        "analysis_status": "frozen spatial-expansion sensitivity pilot",
        "events": public_events,
        "pair_null": pair_summary,
        "control_pair_count": len(null_rows),
        "maximum_measured_vs_fixed_vertex_shift_rsun": max_geometry_shift,
        "measured_pair_order_stable_two_of_three": measured_stable,
        "measured_pair_transport_two_of_three": measured_transport,
        "measured_pair_raw_p_below_0p05": measured_unusual,
        "previous_timing_only_primary_pair_p_raw": previous["primary_pair_p_raw"],
        "previous_timing_only_transport_pass": previous[
            "passed_expanding_transport_compatibility_rule"
        ],
        "verdict": verdict,
        "physical_interpretation": physical_interpretation,
        "no_BH": True,
        "calibration": "Level-0-derived COR1 tB/pB morphology plus pB projected-width proxy",
        "specification": "PDS_20080113_spatially_expanding_ridge_pilot_spec.md",
    }

    map_path = make_map_figure(observed_maps, observed, radii)
    summary_path = make_summary_figure(observed, pair_summary)
    report_path = write_report(summary)
    summary["files"] = {
        "map_figure": map_path.name,
        "summary_figure": summary_path.name,
        "report": report_path.name,
    }
    json_path = OUT / "pds_20080113_spatially_expanding_ridge_results_no_BH.json"
    json_path.write_text(
        json.dumps(summary, indent=2, default=json_default, allow_nan=True),
        encoding="utf-8",
    )
    csv_path = OUT / "pds_20080113_spatially_expanding_ridge_pair_null.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(null_rows[0]))
        writer.writeheader()
        writer.writerows(null_rows)
    payload = {"offsets_min": OFFSETS_MIN}
    for scheme in SCHEMES:
        for candidate_name in CANDIDATES:
            for representation in REPRESENTATIONS:
                row = observed[scheme][candidate_name][representation]
                prefix = f"{scheme}_{candidate_name}_{representation}"
                payload[f"{prefix}_joint_posterior"] = row["joint_posterior"]
                payload[f"{prefix}_tb_curve"] = row["curves"]["total_b"]
                payload[f"{prefix}_pb_curve"] = row["curves"]["pb"]
    npz_path = OUT / "pds_20080113_spatially_expanding_ridge_posteriors.npz"
    np.savez_compressed(npz_path, **payload)

    print(json.dumps({
        "verdict": verdict,
        "pair_primary_p": {
            scheme: pair_summary[scheme]["primary_pair_p_raw"] for scheme in SCHEMES
        },
        "measured_transport_two_of_three": measured_transport,
        "measured_order_two_of_three": measured_stable,
        "max_vertex_shift_rsun": max_geometry_shift,
        "files": [str(report_path), str(json_path), str(csv_path),
                  str(npz_path), str(map_path), str(summary_path)],
    }, indent=2))


if __name__ == "__main__":
    main()
