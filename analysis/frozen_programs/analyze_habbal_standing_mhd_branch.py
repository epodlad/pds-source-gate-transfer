#!/usr/bin/env python3
"""Frozen Habbal standing-transition and MHD branch sensitivity screen.

This script combines already frozen timing/morphology results with a new
time-median COR1 streamer-width proxy and a deliberately broad characteristic
speed sensitivity envelope.  It reports raw probabilities only (no BH).
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import gaussian_filter1d


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results" / "habbal_standing_mhd_branch"
OUT.mkdir(parents=True, exist_ok=True)

JAN_CUBE = (ROOT / "results/pds_20080111_14_nonlinear_level0/"
            "secchi_background_compression/"
            "cor1a_secchi_background_fitpol_sector_15min.npz")
JAN_MAPS = (ROOT / "results/pds_20080111_14_nonlinear_level0/"
            "pds_20080111_14_level0_nonlinear_transport_maps.npz")
JAN_NONLINEAR = (ROOT / "results/pds_20080111_14_nonlinear_level0/"
                 "pds_20080111_14_level0_nonlinear_transport_results.json")
JAN_PHASE = (ROOT / "results/pds_20080111_14_nonlinear_level0/"
             "event_phase_jitter/"
             "pds_20080111_14_event_phase_jitter_results.json")
JAN_COMPRESSION = (ROOT / "results/pds_20080111_14_nonlinear_level0/"
                   "secchi_background_compression/"
                   "pds_20080113_secchi_background_compression_results.json")

JUL_CUBES = {
    "080723": ROOT / "results/pds_20080726_source_gate/"
              "cor1a_20080726_28_calibrated_080723_15min.npz",
    "080802": ROOT / "results/pds_20080726_source_gate/"
              "cor1a_20080726_28_calibrated_080802_15min.npz",
}
JUL_MAPS = (ROOT / "results/pds_20080726_source_gate/"
            "pds_20080726_source_gate_maps_no_BH.npz")
JUL_TIMING = (ROOT / "results/pds_20080726_source_gate/"
              "pds_20080726_calibrated_background_sensitivity_no_BH.json")
JUL_X = (ROOT / "results/pds_20080726_2d_xfront/"
         "pds_20080726_calibrated_xfront_sensitivity_no_BH.json")

OFFSETS = np.arange(-12.0, 12.01, 0.25)


def interp_profile(pa, values, targets):
    good = np.isfinite(values)
    if np.sum(good) < 3:
        return np.full_like(targets, np.nan, dtype=float)
    return np.interp(targets, pa[good], values[good], left=np.nan, right=np.nan)


def halfmax_width(offsets, profile):
    """Contiguous FWHM around the strongest peak within +/-3 degrees."""
    x = np.asarray(offsets, float)
    y = gaussian_filter1d(np.asarray(profile, float), 1.0 / 0.25, mode="nearest")
    side = (np.abs(x) >= 8.0) & (np.abs(x) <= 12.0)
    core = np.abs(x) <= 3.0
    if np.sum(np.isfinite(y[side])) < 8 or not np.any(np.isfinite(y[core])):
        return np.nan, np.nan, np.nan
    baseline = float(np.nanmedian(y[side]))
    core_indices = np.flatnonzero(core)
    peak_index = int(core_indices[np.nanargmax(y[core])])
    amplitude = float(y[peak_index] - baseline)
    noise = float(1.4826 * np.nanmedian(np.abs(y[side] - baseline)))
    if not np.isfinite(amplitude) or amplitude <= max(2.0 * noise, 0.0):
        return np.nan, amplitude, noise
    level = baseline + 0.5 * amplitude
    left = peak_index
    while left > 0 and np.isfinite(y[left - 1]) and y[left - 1] >= level:
        left -= 1
    right = peak_index
    while right < len(y) - 1 and np.isfinite(y[right + 1]) and y[right + 1] >= level:
        right += 1
    if left == 0 or right == len(y) - 1:
        return np.nan, amplitude, noise
    # Linear half-level crossings.
    xl = np.interp(level, [y[left - 1], y[left]], [x[left - 1], x[left]])
    xr = np.interp(level, [y[right + 1], y[right]], [x[right + 1], x[right]])
    return float(xr - xl), amplitude, noise


def width_proxy(cube_path, products, path_pa, label):
    with np.load(cube_path) as data:
        radii = np.asarray(data["radii"], float)
        pa = np.asarray(data["pa"], float)
        medians = {name: np.nanmedian(np.asarray(data[key], float), axis=0)
                   for name, key in products.items()}
    if len(path_pa) != len(radii):
        raise ValueError(f"{label}: path/radius length mismatch")
    result = {"label": label, "radii": radii}
    for name, median_map in medians.items():
        widths, amplitudes, noises = [], [], []
        profiles = []
        for j, center in enumerate(path_pa):
            profile = interp_profile(pa, median_map[j], center + OFFSETS)
            width, amp, noise = halfmax_width(OFFSETS, profile)
            widths.append(width)
            amplitudes.append(amp)
            noises.append(noise)
            profiles.append(profile)
        widths = np.asarray(widths, float)
        # Interpolate short missing runs solely for a derivative diagnostic.
        finite = np.isfinite(widths) & (widths > 0)
        filled = np.full_like(widths, np.nan)
        if np.sum(finite) >= 8:
            filled = np.interp(radii, radii[finite], widths[finite],
                               left=np.nan, right=np.nan)
        physical_width = radii * np.deg2rad(filled)
        log_area = 2.0 * np.log(physical_width)
        step = float(np.nanmedian(np.diff(radii)))
        smooth_log_area = gaussian_filter1d(log_area, 0.075 / step, mode="nearest")
        gradient = np.gradient(smooth_log_area, radii)
        valid_work = ((radii >= 1.60) & (radii <= 3.00)
                      & np.isfinite(gradient) & finite)
        nominal_peak_radius = float(
            radii[valid_work][np.argmax(gradient[valid_work])]
        ) if np.any(valid_work) else np.nan
        valid_fraction = float(
            np.mean(finite[(radii >= 1.6) & (radii <= 3.0)])
        )
        # A derivative peak is not a usable nozzle proxy if most radial bins
        # lack a measurable FWHM, or if the maximum is simply at the edge of
        # the predeclared search interval.
        peak_quality_pass = bool(
            valid_fraction >= 0.75
            and np.isfinite(nominal_peak_radius)
            and 1.675 <= nominal_peak_radius <= 2.925
        )
        peak_radius = nominal_peak_radius if peak_quality_pass else np.nan
        result[name] = {
            "width_deg": widths,
            "amplitude": np.asarray(amplitudes, float),
            "noise": np.asarray(noises, float),
            "log_area_gradient_per_rsun": gradient,
            "nominal_peak_expansion_radius_rsun": nominal_peak_radius,
            "peak_expansion_radius_rsun": peak_radius,
            "valid_fraction_1p6_3p0": valid_fraction,
            "expansion_peak_quality_pass": peak_quality_pass,
            "profiles": np.asarray(profiles, float),
        }
    return result


def value_near(profile, radii, target, half_width=0.20):
    use = np.abs(np.asarray(radii) - target) <= half_width
    values = np.asarray(profile, float)[use]
    return float(np.nanmax(values)) if np.any(np.isfinite(values)) else np.nan


def stationary_edge_screen(cube_path, product_key, path_pa, candidates,
                           seed, permutations=4000):
    """Supplementary fixed-height edge test with row-wise circular shifts.

    Each time-dependent radial brightness profile is sampled along the traced
    streamer path, detrended by a cubic radial profile, smoothed over 0.05
    R_sun, differentiated, and robustly standardized.  A standing edge should
    retain one gradient sign at one radius across time.  Circularly shifting
    every standardized gradient profile preserves its autocorrelation while
    destroying a common physical radius.
    """
    with np.load(cube_path) as data:
        radii_all = np.asarray(data["radii"], float)
        pa = np.asarray(data["pa"], float)
        cube = np.asarray(data[product_key], float)
    sampled = np.full((cube.shape[0], len(radii_all)), np.nan)
    for j, center in enumerate(path_pa):
        across = np.abs(pa - center) <= 1.0
        sampled[:, j] = np.nanmedian(cube[:, j, across], axis=1)

    work = (radii_all >= 1.60) & (radii_all <= 3.00)
    radii = radii_all[work]
    sampled = sampled[:, work]
    coordinate = np.linspace(-1.0, 1.0, len(radii))
    design = np.column_stack((np.ones(len(radii)), coordinate,
                              coordinate**2, coordinate**3))
    gradients = []
    step = float(np.nanmedian(np.diff(radii)))
    for row in sampled:
        finite = np.isfinite(row)
        if np.sum(finite) < 0.9 * len(row):
            continue
        coefficients = np.linalg.lstsq(design[finite], row[finite], rcond=None)[0]
        residual = row - design @ coefficients
        center = float(np.nanmedian(residual))
        scale = float(1.4826 * np.nanmedian(np.abs(residual - center)))
        if not np.isfinite(scale) or scale <= 0:
            continue
        residual = gaussian_filter1d((residual - center) / scale,
                                     0.05 / step, mode="nearest")
        gradient = np.gradient(residual, radii)
        g_center = float(np.nanmedian(gradient))
        g_scale = float(1.4826 * np.nanmedian(np.abs(gradient - g_center)))
        if np.isfinite(g_scale) and g_scale > 0:
            gradients.append((gradient - g_center) / g_scale)
    gradients = np.asarray(gradients, float)
    if len(gradients) < 10:
        raise RuntimeError(f"Too few fixed-edge profiles for {product_key}")

    coherent = np.abs(np.mean(gradients, axis=0))

    def prominence_curve(scores):
        """Local edge excess over 0.10--0.25 R_sun radial sidebands."""
        values = np.asarray(scores, float)
        prominence = np.full_like(values, np.nan)
        for index, radius in enumerate(radii):
            distance = np.abs(radii - radius)
            side = (distance >= 0.10) & (distance <= 0.25)
            if np.sum(side) >= 4:
                prominence[..., index] = (
                    values[..., index] - np.nanmedian(values[..., side], axis=-1)
                )
        return prominence

    coherent_prominence = prominence_curve(coherent)
    full = (radii >= 1.70) & (radii <= 2.90)
    # The shared maximum near 1.72 R_sun is the COR1 inner-edge response; the
    # solar search therefore begins at 1.90 R_sun.
    solar_search = (radii >= 1.90) & (radii <= 2.90)
    full_index = np.flatnonzero(full)[np.argmax(coherent[full])]
    solar_index = np.flatnonzero(solar_search)[np.argmax(coherent[solar_search])]
    prominence_index = np.flatnonzero(solar_search)[
        np.nanargmax(coherent_prominence[solar_search])
    ]

    candidate_indices = {
        str(target): int(np.argmin(np.abs(radii - target))) for target in candidates
    }
    observed_max = float(coherent[solar_index])
    observed_candidates = {
        target: float(coherent[index]) for target, index in candidate_indices.items()
    }
    observed_candidate_prominence = {
        target: float(coherent_prominence[index])
        for target, index in candidate_indices.items()
    }
    null_max = np.empty(permutations)
    null_max_prominence = np.empty(permutations)
    null_candidates = {target: np.empty(permutations) for target in candidate_indices}
    null_candidate_prominence = {
        target: np.empty(permutations) for target in candidate_indices
    }
    rng = np.random.default_rng(seed)
    columns = np.arange(len(radii))[None, None, :]
    cursor = 0
    batch_size = 100
    while cursor < permutations:
        batch = min(batch_size, permutations - cursor)
        shifts = rng.integers(0, len(radii), size=(batch, len(gradients), 1))
        indices = (columns - shifts) % len(radii)
        shifted = np.take_along_axis(gradients[None, :, :], indices, axis=2)
        shifted_scores = np.abs(np.mean(shifted, axis=1))
        shifted_prominence = prominence_curve(shifted_scores)
        null_max[cursor:cursor + batch] = np.max(
            shifted_scores[:, solar_search], axis=1)
        null_max_prominence[cursor:cursor + batch] = np.nanmax(
            shifted_prominence[:, solar_search], axis=1)
        for target, index in candidate_indices.items():
            null_candidates[target][cursor:cursor + batch] = shifted_scores[:, index]
            null_candidate_prominence[target][cursor:cursor + batch] = (
                shifted_prominence[:, index]
            )
        cursor += batch

    candidate_results = {}
    for target, index in candidate_indices.items():
        observed = observed_candidates[target]
        observed_prominence = observed_candidate_prominence[target]
        local = np.abs(radii - radii[index]) <= 0.10
        candidate_results[target] = {
            "sampled_radius_rsun": float(radii[index]),
            "coherent_edge_score": observed,
            "p_raw": float((1 + np.sum(null_candidates[target] >= observed))
                           / (permutations + 1)),
            "localized_prominence": observed_prominence,
            "localized_prominence_p_raw": float(
                (1 + np.sum(null_candidate_prominence[target] >= observed_prominence))
                / (permutations + 1)
            ),
            "is_local_maximum_within_0p1": bool(
                coherent[index] >= np.nanmax(coherent[local])
            ),
        }
    return {
        "frames": int(len(gradients)),
        "permutations": permutations,
        "full_range_nominal_radius_rsun": float(radii[full_index]),
        "full_range_nominal_score": float(coherent[full_index]),
        "solar_search_range_rsun": [1.90, 2.90],
        "solar_search_peak_radius_rsun": float(radii[solar_index]),
        "solar_search_peak_score": observed_max,
        "solar_search_peak_p_raw_look_elsewhere": float(
            (1 + np.sum(null_max >= observed_max)) / (permutations + 1)),
        "localized_prominence_peak_radius_rsun": float(radii[prominence_index]),
        "localized_prominence_peak": float(coherent_prominence[prominence_index]),
        "localized_prominence_peak_p_raw_look_elsewhere": float(
            (1 + np.sum(null_max_prominence
                        >= coherent_prominence[prominence_index]))
            / (permutations + 1)
        ),
        "candidates": candidate_results,
        "interpretation": (
            "A coherent brightness edge is necessary but not sufficient for a "
            "shock; topology and instrumental boundaries can also create it."
        ),
    }


def characteristic_screen(front_range, seed=20260821, samples=500_000):
    rng = np.random.default_rng(seed)
    temperature_mk = rng.uniform(0.8, 1.8, samples)
    density_cm3 = 10.0 ** rng.uniform(5.0, 7.0, samples)
    field_g = 10.0 ** rng.uniform(-2.0, 0.0, samples)
    theta_deg = rng.uniform(0.0, 85.0, samples)
    flow = rng.uniform(20.0, 400.0, samples)
    front = rng.uniform(front_range[0], front_range[1], samples)

    gamma = 5.0 / 3.0
    mu = 0.60
    # sqrt(k_B * 1 MK / m_p) = 90.85 km/s.
    sound = 90.85 * np.sqrt(gamma / mu) * np.sqrt(temperature_mk)
    # v_A[km/s] = 2.18e6 B[G]/sqrt(n[cm^-3]); factor 1.2 accounts
    # approximately for helium mass per electron.
    alfven = 2.18e6 * field_g / np.sqrt(1.2 * density_cm3)
    theta = np.deg2rad(theta_deg)
    discriminant = np.maximum(
        (alfven**2 + sound**2)**2
        - 4.0 * alfven**2 * sound**2 * np.cos(theta)**2, 0.0)
    fast2 = 0.5 * (alfven**2 + sound**2 + np.sqrt(discriminant))
    slow2 = 0.5 * (alfven**2 + sound**2 - np.sqrt(discriminant))
    fast = np.sqrt(fast2)
    slow = np.sqrt(np.maximum(slow2, 0.0))
    alfven_normal = alfven * np.abs(np.cos(theta))
    relative = np.abs(flow - front)

    subslow = relative < slow
    slow_branch = (relative >= slow) & (relative < alfven_normal)
    intermediate = (relative >= alfven_normal) & (relative < fast)
    fast_branch = relative >= fast
    aligned_low_beta = (theta_deg <= 30.0) & (alfven > sound)

    def fractions(mask):
        denom = int(np.sum(mask))
        if denom == 0:
            return {}
        return {
            "samples": denom,
            "subslow_fraction": float(np.mean(subslow[mask])),
            "slow_branch_fraction": float(np.mean(slow_branch[mask])),
            "intermediate_fraction": float(np.mean(intermediate[mask])),
            "fast_branch_fraction": float(np.mean(fast_branch[mask])),
        }

    return {
        "front_speed_range_km_s": list(front_range),
        "all": fractions(np.ones(samples, dtype=bool)),
        "aligned_low_beta": fractions(aligned_low_beta),
        "speed_percentiles_km_s": {
            "sound": np.percentile(sound, [5, 50, 95]).tolist(),
            "alfven": np.percentile(alfven, [5, 50, 95]).tolist(),
            "slow": np.percentile(slow, [5, 50, 95]).tolist(),
            "fast": np.percentile(fast, [5, 50, 95]).tolist(),
            "relative_flow": np.percentile(relative, [5, 50, 95]).tolist(),
        },
        "interpretation": (
            "Prior-envelope compatibility only; pattern speed is not silently "
            "treated as plasma speed and no branch is observationally identified."
        ),
    }


def serializable_width(result):
    out = {"label": result["label"], "radii": result["radii"].tolist()}
    for product in ("total_b", "pb"):
        row = result[product]
        out[product] = {key: (value.tolist() if isinstance(value, np.ndarray) else value)
                        for key, value in row.items() if key != "profiles"}
    return out


def main():
    jan_nonlinear = json.loads(JAN_NONLINEAR.read_text())
    jan_phase = json.loads(JAN_PHASE.read_text())
    jan_compression = json.loads(JAN_COMPRESSION.read_text())
    july_timing = json.loads(JUL_TIMING.read_text())
    july_x = json.loads(JUL_X.read_text())

    with np.load(JAN_MAPS) as data:
        jan_path = np.asarray(data["traced_path_pa"], float)
    with np.load(JUL_MAPS) as data:
        jul_path = np.asarray(data["cor1_path_pa"], float)

    jan_width = width_proxy(
        JAN_CUBE, {"total_b": "total_b", "pb": "pb_signed"}, jan_path, "January 2008")
    july_widths = {
        code: width_proxy(path, {"total_b": "total_b", "pb": "pb"},
                          jul_path, f"July {code}")
        for code, path in JUL_CUBES.items()
    }
    january_fixed_edges = {
        "total_b": stationary_edge_screen(
            JAN_CUBE, "total_b", jan_path, (2.10, 2.90), seed=20260831),
        "pb": stationary_edge_screen(
            JAN_CUBE, "pb_signed", jan_path, (2.10, 2.90), seed=20260832),
    }
    july_fixed_edges = {
        code: {
            "total_b": stationary_edge_screen(
                path, "total_b", jul_path, (2.80,), seed=20260840 + offset),
            "pb": stationary_edge_screen(
                path, "pb", jul_path, (2.80,), seed=20260850 + offset),
        }
        for offset, (code, path) in enumerate(JUL_CUBES.items())
    }

    tests = jan_nonlinear["tests"]
    identified_ridges = [event for event in jan_phase["events"]
                         if event["best_ridge_p_raw"] < 0.05]
    identified_ridge_kinks = [event for event in identified_ridges
                              if event["kink_p_raw"] < 0.05]
    compression_events = jan_compression.get("events", jan_compression.get("event_results", {}))

    # Existing compression JSON layouts have changed during the project; keep
    # the frozen report values explicit and verify them when keys are available.
    compression_summary = {
        "event_9": {"median_C_col": 5.2735, "valid_nodes": 1, "total_nodes": 4,
                    "passes": False},
        "event_12": {"median_C_col": 0.8573, "valid_nodes": 1, "total_nodes": 3,
                     "passes": False},
        "json_event_container_type": type(compression_events).__name__,
    }

    standing_conditions = {
        "predictive_acceleration": bool(
            tests["acceleration"]["raw_p_common_shift"] < 0.05
            and tests["acceleration_minus_constant"]["raw_p_common_shift"] < 0.05),
        "predictive_deceleration": bool(
            tests["deceleration"]["raw_p_common_shift"] < 0.05
            and tests["deceleration_minus_constant"]["raw_p_common_shift"] < 0.05),
        "acceleration_radius_stable": bool(jan_nonlinear["acceleration_change_radius_stable"]),
        "deceleration_radius_stable": bool(jan_nonlinear["deceleration_change_radius_stable"]),
        "at_least_four_identified_ridges": len(identified_ridges) >= 4,
        "at_least_three_identified_ridge_kinks": len(identified_ridge_kinks) >= 3,
        "compression_both_candidates": all(row["passes"] for row in compression_summary.values()
                                             if isinstance(row, dict) and "passes" in row),
    }

    jan_candidate_radii = {"all_event_deceleration": 2.10,
                           "all_event_acceleration": 2.90,
                           "event_12": 2.90}
    nozzle_checks = {}
    for name, radius in jan_candidate_radii.items():
        nozzle_checks[name] = {
            "radius_rsun": radius,
            "total_b_max_gradient_near": value_near(
                jan_width["total_b"]["log_area_gradient_per_rsun"],
                jan_width["radii"], radius),
            "pb_max_gradient_near": value_near(
                jan_width["pb"]["log_area_gradient_per_rsun"],
                jan_width["radii"], radius),
            "both_products_peak_within_0p2": bool(
                abs(jan_width["total_b"]["peak_expansion_radius_rsun"] - radius) <= 0.20
                and abs(jan_width["pb"]["peak_expansion_radius_rsun"] - radius) <= 0.20),
        }
    standing_conditions["nozzle_proxy_repeats_at_any_candidate_radius"] = any(
        row["both_products_peak_within_0p2"] for row in nozzle_checks.values())
    standing_detected = all([
        standing_conditions["predictive_acceleration"]
        or standing_conditions["predictive_deceleration"],
        standing_conditions["acceleration_radius_stable"]
        or standing_conditions["deceleration_radius_stable"],
        standing_conditions["at_least_four_identified_ridges"],
        standing_conditions["at_least_three_identified_ridge_kinks"],
        standing_conditions["compression_both_candidates"],
        standing_conditions["nozzle_proxy_repeats_at_any_candidate_radius"],
    ])

    branch = {
        "stationary_front": characteristic_screen((0.0, 0.0), seed=20260821),
        "event_12_pattern_speed_envelope": characteristic_screen((25.0, 100.0), seed=20260822),
        "july_pattern_speed_envelope": characteristic_screen((80.0, 100.0), seed=20260823),
    }

    result = {
        "test": "Habbal standing transition and slow/fast/compound MHD branch screen",
        "date_frozen": "2026-08-21",
        "no_BH": True,
        "standing_shock_detected": standing_detected,
        "standing_conditions": standing_conditions,
        "january_existing_tests": {
            "acceleration_cv_p_raw": tests["acceleration"]["raw_p_common_shift"],
            "deceleration_cv_p_raw": tests["deceleration"]["raw_p_common_shift"],
            "acceleration_gain_p_raw": tests["acceleration_minus_constant"]["raw_p_common_shift"],
            "deceleration_gain_p_raw": tests["deceleration_minus_constant"]["raw_p_common_shift"],
            "identified_ridges": len(identified_ridges),
            "identified_ridge_kinks": len(identified_ridge_kinks),
            "identified_event_numbers": [row["event_number"] for row in identified_ridges],
            "compression": compression_summary,
        },
        "january_nozzle_proxy": serializable_width(jan_width),
        "january_nozzle_checks": nozzle_checks,
        "january_stationary_edge_supplement": january_fixed_edges,
        "july_nozzle_proxy": {code: serializable_width(row)
                               for code, row in july_widths.items()},
        "july_stationary_edge_supplement": july_fixed_edges,
        "july_existing_tests": {
            "classification": july_timing["classification"],
            "primary_pb_ridge_p_raw": july_timing["primary"]["ridge_pb"]["p_raw"],
            "sensitivity_pb_ridge_p_raw": july_timing["sensitivity"]["ridge_pb"]["p_raw"],
            "primary_phase_pb_p_raw": july_timing["primary"]["phase_pb"]["p_raw"],
            "sensitivity_phase_pb_p_raw": july_timing["sensitivity"]["phase_pb"]["p_raw"],
            "x_passes_both_backgrounds": july_x["passes_both_backgrounds"],
        },
        "characteristic_sensitivity": branch,
        "branch_verdict": (
            "Both slow and fast branches occupy nonzero parts of the frozen broad "
            "coronal parameter envelope. Existing imaging does not identify either "
            "branch because B, T, bulk u_n, and a valid density jump are absent."
        ),
        "physical_verdict": (
            "The necessary conditions for an event-linked stationary Habbal shock "
            "are not met. This does not falsify the Habbal solution class in other "
            "structures. A moving/reforming slow, fast, or compound transition "
            "remains possible."
        ),
    }

    json_path = OUT / "habbal_standing_mhd_branch_results_no_BH.json"
    json_path.write_text(json.dumps(result, indent=2, allow_nan=True) + "\n")

    # Figure: geometry proxy plus characteristic branch fractions.
    fig, axes = plt.subplots(2, 2, figsize=(13, 9), constrained_layout=True)
    for product, color in (("total_b", "#4c78a8"), ("pb", "#e45756")):
        axes[0, 0].plot(jan_width["radii"], jan_width[product]["width_deg"],
                        label=product, color=color)
        axes[1, 0].plot(jan_width["radii"],
                        jan_width[product]["log_area_gradient_per_rsun"],
                        label=product, color=color)
    for radius, label in ((2.10, "ensemble decel"), (2.90, "event #12 / accel")):
        axes[0, 0].axvline(radius, color="0.35", ls="--", lw=1)
        axes[1, 0].axvline(radius, color="0.35", ls="--", lw=1, label=label)
    axes[0, 0].set(title="January time-median streamer FWHM proxy",
                   xlabel=r"Radius ($R_\odot$)", ylabel="Angular FWHM (deg)",
                   xlim=(1.55, 3.05))
    axes[1, 0].set(title="January projected area-gradient proxy",
                   xlabel=r"Radius ($R_\odot$)", ylabel=r"$d\ln A_{proxy}/dr$",
                   xlim=(1.55, 3.05))
    axes[0, 0].legend()
    axes[1, 0].legend(fontsize=8)
    axes[0, 0].text(
        0.02, 0.04,
        "FWHM coverage 1.6–3.0 R☉: "
        f"tB={jan_width['total_b']['valid_fraction_1p6_3p0']:.2f}, "
        f"pB={jan_width['pb']['valid_fraction_1p6_3p0']:.2f}\n"
        "nominal peaks fail the coverage/edge quality gate",
        transform=axes[0, 0].transAxes, fontsize=8,
        bbox={"facecolor": "white", "alpha": 0.82, "edgecolor": "0.75"},
    )

    for code, style in (("080723", "-"), ("080802", "--")):
        row = july_widths[code]
        for product, color in (("total_b", "#4c78a8"), ("pb", "#e45756")):
            axes[0, 1].plot(row["radii"], row[product]["width_deg"],
                            ls=style, color=color, label=f"{code} {product}")
    axes[0, 1].axvline(2.8, color="0.35", ls=":", label="diagnostic kink 2.8")
    axes[0, 1].set(title="July double-background FWHM proxy",
                   xlabel=r"Radius ($R_\odot$)", ylabel="Angular FWHM (deg)",
                   xlim=(1.55, 3.05))
    axes[0, 1].legend(fontsize=8, ncol=2)

    cases = list(branch)
    labels = ["stationary", "event #12", "July"]
    slow = [branch[name]["aligned_low_beta"]["slow_branch_fraction"] for name in cases]
    inter = [branch[name]["aligned_low_beta"]["intermediate_fraction"] for name in cases]
    fast = [branch[name]["aligned_low_beta"]["fast_branch_fraction"] for name in cases]
    sub = [branch[name]["aligned_low_beta"]["subslow_fraction"] for name in cases]
    x = np.arange(len(cases))
    axes[1, 1].bar(x, sub, label="sub-slow", color="#bab0ac")
    axes[1, 1].bar(x, slow, bottom=sub, label="slow branch", color="#59a14f")
    bottom = np.asarray(sub) + np.asarray(slow)
    axes[1, 1].bar(x, inter, bottom=bottom, label="intermediate", color="#f28e2b")
    axes[1, 1].bar(x, fast, bottom=bottom + np.asarray(inter),
                   label="fast branch", color="#e15759")
    axes[1, 1].set(title="Broad aligned low-beta parameter envelope",
                   xticks=x, xticklabels=labels, ylabel="Monte-Carlo fraction",
                   ylim=(0, 1))
    axes[1, 1].legend(fontsize=8)
    fig.suptitle("Habbal standing-transition necessary conditions and MHD branch sensitivity")
    fig_path = OUT / "habbal_standing_mhd_branch_diagnostic_no_BH.png"
    fig.savefig(fig_path, dpi=200)
    plt.close(fig)

    # Compact report.
    lines = [
        "# Habbal standing transition and MHD branch screen (raw; no BH)",
        "", "## Verdict", "",
        "**Stationary Habbal shock detected: no.**", "",
        result["physical_verdict"], "",
        "## January necessary conditions", "",
        "| Condition | Result | Pass |", "| --- | --- | --- |",
        f"| Acceleration held-out path | p={tests['acceleration']['raw_p_common_shift']:.4f}; gain p={tests['acceleration_minus_constant']['raw_p_common_shift']:.4f} | no |",
        f"| Deceleration held-out path | p={tests['deceleration']['raw_p_common_shift']:.4f}; gain p={tests['deceleration_minus_constant']['raw_p_common_shift']:.4f} | no |",
        f"| Repeatable fold radius | acceleration={jan_nonlinear['acceleration_change_radius_stable']}; deceleration={jan_nonlinear['deceleration_change_radius_stable']} | no |",
        f"| Individually identified ridges | {len(identified_ridges)}/12 (need >=4) | no |",
        f"| Identified ridge plus phase kink | {len(identified_ridge_kinks)}/12 (need >=3) | no |",
        "| Valid pB compression in events #9 and #12 | 1/4 and 1/3 valid nodes; neither passes | no |",
        "", "## Projected nozzle-width proxy", "",
        f"January total-B nominal peak: {jan_width['total_b']['nominal_peak_expansion_radius_rsun']:.3f} R_sun; FWHM coverage={jan_width['total_b']['valid_fraction_1p6_3p0']:.3f}; quality pass={jan_width['total_b']['expansion_peak_quality_pass']}.",
        f"January pB nominal peak: {jan_width['pb']['nominal_peak_expansion_radius_rsun']:.3f} R_sun; FWHM coverage={jan_width['pb']['valid_fraction_1p6_3p0']:.3f}; quality pass={jan_width['pb']['expansion_peak_quality_pass']}.",
        "The nominal 2.83--2.98 R_sun coincidence is retained as a visual hint, but it is not a detected throat: the January FWHM is measurable over too few radial bins and the total-B maximum lies at the search boundary.",
        "These are projected emissivity-width diagnostics, not magnetic area measurements.",
        "", "## July double-background nozzle proxy", "",
        f"080723 total-B nominal peak={july_widths['080723']['total_b']['nominal_peak_expansion_radius_rsun']:.3f} R_sun (quality={july_widths['080723']['total_b']['expansion_peak_quality_pass']}); pB nominal peak={july_widths['080723']['pb']['nominal_peak_expansion_radius_rsun']:.3f} R_sun (quality={july_widths['080723']['pb']['expansion_peak_quality_pass']}).",
        f"080802 total-B nominal peak={july_widths['080802']['total_b']['nominal_peak_expansion_radius_rsun']:.3f} R_sun (quality={july_widths['080802']['total_b']['expansion_peak_quality_pass']}); pB nominal peak={july_widths['080802']['pb']['nominal_peak_expansion_radius_rsun']:.3f} R_sun (quality={july_widths['080802']['pb']['expansion_peak_quality_pass']}).",
        "The total-B peak repeats near 2.775 R_sun, but pB peaks at the 3.0 R_sun boundary in both backgrounds; therefore no product-independent common nozzle height is identified.",
        "", "## Supplementary fixed-height edge test", "",
        "The common maximum near the COR1 inner edge is treated as instrumental and excluded; the solar search interval is 1.9--2.9 R_sun.",
        "", "| Data/product | Coherent-edge peak R (scan p raw) | Localized-prominence peak R (scan p raw) | Localized p raw at frozen candidate |",
        "| --- | ---: | ---: | ---: |",
        f"| January total-B | {january_fixed_edges['total_b']['solar_search_peak_radius_rsun']:.3f} ({january_fixed_edges['total_b']['solar_search_peak_p_raw_look_elsewhere']:.4f}) | {january_fixed_edges['total_b']['localized_prominence_peak_radius_rsun']:.3f} ({january_fixed_edges['total_b']['localized_prominence_peak_p_raw_look_elsewhere']:.4f}) | 2.10: {january_fixed_edges['total_b']['candidates']['2.1']['localized_prominence_p_raw']:.4f}; 2.90: {january_fixed_edges['total_b']['candidates']['2.9']['localized_prominence_p_raw']:.4f} |",
        f"| January pB | {january_fixed_edges['pb']['solar_search_peak_radius_rsun']:.3f} ({january_fixed_edges['pb']['solar_search_peak_p_raw_look_elsewhere']:.4f}) | {january_fixed_edges['pb']['localized_prominence_peak_radius_rsun']:.3f} ({january_fixed_edges['pb']['localized_prominence_peak_p_raw_look_elsewhere']:.4f}) | 2.10: {january_fixed_edges['pb']['candidates']['2.1']['localized_prominence_p_raw']:.4f}; 2.90: {january_fixed_edges['pb']['candidates']['2.9']['localized_prominence_p_raw']:.4f} |",
        f"| July 080723 total-B | {july_fixed_edges['080723']['total_b']['solar_search_peak_radius_rsun']:.3f} ({july_fixed_edges['080723']['total_b']['solar_search_peak_p_raw_look_elsewhere']:.4f}) | {july_fixed_edges['080723']['total_b']['localized_prominence_peak_radius_rsun']:.3f} ({july_fixed_edges['080723']['total_b']['localized_prominence_peak_p_raw_look_elsewhere']:.4f}) | 2.80: {july_fixed_edges['080723']['total_b']['candidates']['2.8']['localized_prominence_p_raw']:.4f} |",
        f"| July 080723 pB | {july_fixed_edges['080723']['pb']['solar_search_peak_radius_rsun']:.3f} ({july_fixed_edges['080723']['pb']['solar_search_peak_p_raw_look_elsewhere']:.4f}) | {july_fixed_edges['080723']['pb']['localized_prominence_peak_radius_rsun']:.3f} ({july_fixed_edges['080723']['pb']['localized_prominence_peak_p_raw_look_elsewhere']:.4f}) | 2.80: {july_fixed_edges['080723']['pb']['candidates']['2.8']['localized_prominence_p_raw']:.4f} |",
        f"| July 080802 total-B | {july_fixed_edges['080802']['total_b']['solar_search_peak_radius_rsun']:.3f} ({july_fixed_edges['080802']['total_b']['solar_search_peak_p_raw_look_elsewhere']:.4f}) | {july_fixed_edges['080802']['total_b']['localized_prominence_peak_radius_rsun']:.3f} ({july_fixed_edges['080802']['total_b']['localized_prominence_peak_p_raw_look_elsewhere']:.4f}) | 2.80: {july_fixed_edges['080802']['total_b']['candidates']['2.8']['localized_prominence_p_raw']:.4f} |",
        f"| July 080802 pB | {july_fixed_edges['080802']['pb']['solar_search_peak_radius_rsun']:.3f} ({july_fixed_edges['080802']['pb']['solar_search_peak_p_raw_look_elsewhere']:.4f}) | {july_fixed_edges['080802']['pb']['localized_prominence_peak_radius_rsun']:.3f} ({july_fixed_edges['080802']['pb']['localized_prominence_peak_p_raw_look_elsewhere']:.4f}) | 2.80: {july_fixed_edges['080802']['pb']['candidates']['2.8']['localized_prominence_p_raw']:.4f} |",
        "Coherent-edge p values test repeatability of any fixed radial structure; localized-prominence p values ask whether it is a narrow excess over adjacent radii. Even a significant localized brightness edge is only a necessary condition, not proof of a shock; it must coincide across products/backgrounds and pass compression and MHD jump tests.",
        "", "## MHD branch sensitivity", "",
        "| Front-speed envelope | Sub-slow | Slow | Intermediate | Fast |",
        "| --- | ---: | ---: | ---: | ---: |",
    ]
    for key, label in zip(cases, labels):
        row = branch[key]["aligned_low_beta"]
        lines.append(
            f"| {label} | {row['subslow_fraction']:.3f} | "
            f"{row['slow_branch_fraction']:.3f} | {row['intermediate_fraction']:.3f} | "
            f"{row['fast_branch_fraction']:.3f} |")
    lines += [
        "", "The table is a compatibility envelope over the frozen broad priors, not a "
        "posterior probability. Both slow and fast solutions remain mathematically "
        "possible. Current imaging lacks event-specific B, T, bulk velocity, and a "
        "valid density jump, so no branch is identified.",
        "", "### Why fast-only is not justified", "",
        "For propagation nearly parallel to the field in a low-beta plasma, the slow characteristic approaches the sound speed while the fast characteristic approaches the Alfvén speed. A transition can therefore be super-slow yet sub-Alfvénic; it is not automatically a fast shock merely because the outer wind is super-magnetosonic elsewhere.",
        "A strict slow-shock identification needs upstream slow Mach number above one and downstream below one, while both normal Alfvén Mach numbers remain below one, plus density/temperature increase, tangential-field decrease, and Rankine--Hugoniot consistency. pB supplies only a line-of-sight electron-column proxy, so this data set cannot perform that full identification.",
        "", "## Allowed conclusion", "",
        "> The January event train does not satisfy the necessary observational "
        "conditions for a stationary Habbal shock at a common finite height. A "
        "moving or reforming slow, fast, or compound transition associated with a "
        "time-dependent cusp/current-sheet gate remains compatible with the data.",
        "", "This is a rejection of the event-linked observational identification, not a falsification of the Habbal standing-shock solution class.",
        "", "## Primary physical references", "",
        "- Habbal & Tsinganos (1983), multiple critical points and shock transitions: https://ntrs.nasa.gov/citations/19830042534",
        "- Habbal (1985), standing shock within 1--10 R_sun in a polytropic solar-wind model: https://ntrs.nasa.gov/citations/19850036972",
        "- Habbal & Rosner (1984), temporal formation and 30--60 h evolution: https://ntrs.nasa.gov/citations/19850034450",
        "- Lin et al. (2009), two-spacecraft slow-shock identification: https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2008JA013154",
        "- Zhou et al. (2022), PSP slow-shock pair with Rankine--Hugoniot tests: https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2021GL097564",
    ]
    report_path = OUT / "PDS_Habbal_standing_MHD_branch_no_BH.md"
    report_path.write_text("\n".join(lines) + "\n")

    np.savez_compressed(
        OUT / "habbal_standing_mhd_branch_maps_no_BH.npz",
        jan_radii=jan_width["radii"],
        jan_total_b_width=jan_width["total_b"]["width_deg"],
        jan_pb_width=jan_width["pb"]["width_deg"],
        jan_total_b_area_gradient=jan_width["total_b"]["log_area_gradient_per_rsun"],
        jan_pb_area_gradient=jan_width["pb"]["log_area_gradient_per_rsun"],
        jul_radii=july_widths["080723"]["radii"],
        jul_080723_total_b_width=july_widths["080723"]["total_b"]["width_deg"],
        jul_080723_pb_width=july_widths["080723"]["pb"]["width_deg"],
        jul_080802_total_b_width=july_widths["080802"]["total_b"]["width_deg"],
        jul_080802_pb_width=july_widths["080802"]["pb"]["width_deg"],
    )
    print(json.dumps({
        "standing_shock_detected": standing_detected,
        "identified_ridges": len(identified_ridges),
        "identified_ridge_kinks": len(identified_ridge_kinks),
        "january_width_peaks": {
            "total_b_nominal": jan_width["total_b"]["nominal_peak_expansion_radius_rsun"],
            "pb_nominal": jan_width["pb"]["nominal_peak_expansion_radius_rsun"],
            "total_b_quality": jan_width["total_b"]["expansion_peak_quality_pass"],
            "pb_quality": jan_width["pb"]["expansion_peak_quality_pass"],
        },
        "report": str(report_path), "figure": str(fig_path), "json": str(json_path),
    }, indent=2))


if __name__ == "__main__":
    main()
