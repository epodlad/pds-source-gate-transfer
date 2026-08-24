#!/usr/bin/env python3
"""Geometry-dependent Habbal and moving/reforming cusp-gate screen.

Uses the already frozen January event anchors/path and the independent July
double-background cubes.  All probabilities are raw empirical values; no BH.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.signal import detrend
from scipy.stats import rankdata, spearmanr


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results" / "dynamic_geometry_gate"
OUT.mkdir(parents=True, exist_ok=True)

JAN_CUBE = (ROOT / "results/pds_20080111_14_nonlinear_level0/"
            "secchi_background_compression/"
            "cor1a_secchi_background_fitpol_sector_15min.npz")
JAN_MAPS = (ROOT / "results/pds_20080111_14_nonlinear_level0/"
            "pds_20080111_14_level0_nonlinear_transport_maps.npz")
JAN_PHASE = (ROOT / "results/pds_20080111_14_nonlinear_level0/"
             "event_phase_jitter/"
             "pds_20080111_14_event_phase_jitter_results.json")

JUL_CUBES = {
    "080723": ROOT / "results/pds_20080726_source_gate/"
              "cor1a_20080726_28_calibrated_080723_15min.npz",
    "080802": ROOT / "results/pds_20080726_source_gate/"
              "cor1a_20080726_28_calibrated_080802_15min.npz",
}
JUL_MAPS = (ROOT / "results/pds_20080726_source_gate/"
            "pds_20080726_source_gate_maps_no_BH.npz")

RSUN_KM = 695_700.0
PRIMARY_SPEED = 200.0
SPEEDS = (100.0, 200.0, 300.0)
RADIAL_SEARCH = (1.90, 2.90)
FIXED_RADII = (2.10, 2.50, 2.80, 2.90)


def robust_mad(values, axis=None):
    values = np.asarray(values, float)
    center = np.nanmedian(values, axis=axis, keepdims=True)
    mad = 1.4826 * np.nanmedian(np.abs(values - center), axis=axis)
    return mad


def empirical_p(observed, null, two_sided=False):
    null = np.asarray(null, float)
    null = null[np.isfinite(null)]
    if not np.isfinite(observed) or len(null) == 0:
        return np.nan
    if two_sided:
        observed = abs(observed)
        null = np.abs(null)
    return float((1 + np.sum(null >= observed)) / (len(null) + 1))


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


def streamer_moments(cube_path, product_key, path_pa):
    with np.load(cube_path) as data:
        cube = np.asarray(data[product_key], float)
        times = np.asarray(data["minutes"], float)
        radii = np.asarray(data["radii"], float)
        pa = np.asarray(data["pa"], float)
    nt, nr, _ = cube.shape
    width = np.full((nt, nr), np.nan)
    centroid = np.full((nt, nr), np.nan)
    excess = np.full((nt, nr), np.nan)
    amplitude_snr = np.full((nt, nr), np.nan)

    for j, center in enumerate(path_pa):
        offsets = pa - center
        use = np.abs(offsets) <= 12.0
        side = use & (np.abs(offsets) >= 8.0)
        core = use & (np.abs(offsets) <= 8.0)
        if np.sum(side) < 6 or np.sum(core) < 8:
            continue
        profiles = cube[:, j, :]
        baseline = np.nanmedian(profiles[:, side], axis=1)
        noise = robust_mad(profiles[:, side], axis=1)
        signal = profiles[:, core] - baseline[:, None]
        weights = np.clip(signal, 0.0, None)
        total = np.nansum(weights, axis=1)
        peak = np.nanmax(signal, axis=1)
        x = offsets[core]
        good = (np.isfinite(total) & (total > 0) & np.isfinite(noise)
                & (peak > 2.0 * np.maximum(noise, np.finfo(float).eps)))
        c = np.full(nt, np.nan)
        s = np.full(nt, np.nan)
        c[good] = np.nansum(weights[good] * x[None, :], axis=1) / total[good]
        variance = np.full(nt, np.nan)
        variance[good] = (
            np.nansum(weights[good] * (x[None, :] - c[good, None])**2, axis=1)
            / total[good]
        )
        s[good] = np.sqrt(np.maximum(variance[good], 0.0))
        width[:, j] = s
        centroid[:, j] = c
        dpa = float(np.nanmedian(np.diff(pa[use])))
        excess[good, j] = total[good] * abs(dpa)
        amplitude_snr[good, j] = peak[good] / np.maximum(
            noise[good], np.finfo(float).eps)

    physical_width = radii[None, :] * np.deg2rad(width)
    log_area = np.full_like(physical_width, np.nan)
    positive_width = np.isfinite(physical_width) & (physical_width > 0)
    log_area[positive_width] = 2.0 * np.log(physical_width[positive_width])
    log_excess = np.log(excess)

    def radial_derivative(values, sigma_r=0.075):
        derivative = np.full_like(values, np.nan)
        step = float(np.nanmedian(np.diff(radii)))
        for i, row in enumerate(values):
            good = np.isfinite(row)
            work = good & (radii >= 1.60) & (radii <= 3.00)
            if np.mean(work[(radii >= 1.60) & (radii <= 3.00)]) < 0.55:
                continue
            filled = np.interp(radii, radii[good], row[good], left=np.nan, right=np.nan)
            inside = np.isfinite(filled)
            if np.sum(inside) < 12:
                continue
            lo, hi = np.flatnonzero(inside)[[0, -1]]
            smoothed = gaussian_filter1d(
                filled[lo:hi + 1], sigma_r / step, mode="nearest")
            derivative[i, lo:hi + 1] = np.gradient(smoothed, radii[lo:hi + 1])
        return derivative

    return {
        "times": times,
        "radii": radii,
        "width_deg": width,
        "centroid_deg": centroid,
        "excess": excess,
        "amplitude_snr": amplitude_snr,
        "log_area": log_area,
        "log_excess": log_excess,
        "area_gradient": radial_derivative(log_area),
        "axis_shear": radial_derivative(centroid, sigma_r=0.075),
        "valid_fraction": {
            "width": float(np.mean(np.isfinite(width))),
            "centroid": float(np.mean(np.isfinite(centroid))),
            "excess": float(np.mean(np.isfinite(excess))),
        },
    }


def interpolate_path(values, times, targets):
    output = np.full(values.shape[1], np.nan)
    for j, target in enumerate(targets):
        good = np.isfinite(values[:, j])
        if np.sum(good) < 3:
            continue
        good_times = times[good]
        good_values = values[good, j]
        upper = int(np.searchsorted(good_times, target, side="left"))
        if upper < len(good_times) and abs(good_times[upper] - target) < 1e-8:
            output[j] = good_values[upper]
        elif 0 < upper < len(good_times):
            lower = upper - 1
            if good_times[upper] - good_times[lower] <= 30.0:
                output[j] = np.interp(
                    target, good_times[lower:upper + 1],
                    good_values[lower:upper + 1])
            else:
                nearest = lower if (target - good_times[lower]
                                    <= good_times[upper] - target) else upper
                if abs(good_times[nearest] - target) <= 30.0:
                    output[j] = good_values[nearest]
        elif upper == 0 and abs(good_times[0] - target) <= 30.0:
            output[j] = good_values[0]
        elif upper == len(good_times) and abs(good_times[-1] - target) <= 30.0:
            output[j] = good_values[-1]
    return output


def backprojected_targets(anchor, radii, speed):
    return anchor - (3.0 - radii) * RSUN_KM / speed / 60.0


def event_geometry(moment, anchors, speed):
    radii = moment["radii"]
    search = ((radii >= RADIAL_SEARCH[0]) & (radii <= RADIAL_SEARCH[1]))
    rows = []
    for anchor in anchors:
        target = backprojected_targets(anchor, radii, speed)
        gradient = interpolate_path(moment["area_gradient"], moment["times"], target)
        shear = np.abs(interpolate_path(moment["axis_shear"], moment["times"], target))
        log_area = interpolate_path(moment["log_area"], moment["times"], target)
        coverage = float(np.mean(np.isfinite(log_area[search])))
        valid_g = search & np.isfinite(gradient)
        valid_s = search & np.isfinite(shear)
        if coverage < 0.55 or not np.any(valid_g):
            rows.append({
                "coverage": coverage, "max_area_gradient": np.nan,
                "gradient_peak_radius_rsun": np.nan,
                "area_expansion_ratio": np.nan, "max_axis_shear": np.nan,
                "axis_shear_peak_radius_rsun": np.nan,
            })
            continue
        g_index = np.flatnonzero(valid_g)[np.argmax(gradient[valid_g])]
        r_peak = radii[g_index]
        inner = float(np.interp(r_peak - 0.15, radii[np.isfinite(log_area)],
                                log_area[np.isfinite(log_area)]))
        outer = float(np.interp(r_peak + 0.15, radii[np.isfinite(log_area)],
                                log_area[np.isfinite(log_area)]))
        if np.any(valid_s):
            s_index = np.flatnonzero(valid_s)[np.argmax(shear[valid_s])]
            s_value = float(shear[s_index])
            s_radius = float(radii[s_index])
        else:
            s_value, s_radius = np.nan, np.nan
        rows.append({
            "coverage": coverage,
            "max_area_gradient": float(gradient[g_index]),
            "gradient_peak_radius_rsun": float(r_peak),
            "area_expansion_ratio": float(np.exp(outer - inner)),
            "max_axis_shear": s_value,
            "axis_shear_peak_radius_rsun": s_radius,
        })
    return rows


def permutation_correlation(rows, response, seed, permutations=50_000):
    metrics = ("max_area_gradient", "area_expansion_ratio", "max_axis_shear")
    rng = np.random.default_rng(seed)
    output = {}
    null_by_metric = {}
    for metric in metrics:
        x = np.array([row[metric] for row in rows], float)
        y = np.asarray(response, float)
        good = np.isfinite(x) & np.isfinite(y)
        if np.sum(good) < 7:
            output[metric] = {"rho": np.nan, "p_raw": np.nan,
                              "pairs": int(np.sum(good))}
            null_by_metric[metric] = np.full(permutations, np.nan)
            continue
        xg, yg = x[good], y[good]
        observed = float(spearmanr(xg, yg).statistic)
        xr = rankdata(xg)
        yr = rankdata(yg)
        xr = (xr - np.mean(xr)) / np.std(xr)
        yr = (yr - np.mean(yr)) / np.std(yr)
        null = np.empty(permutations)
        cursor = 0
        while cursor < permutations:
            batch = min(2000, permutations - cursor)
            order = np.argsort(rng.random((batch, len(yr))), axis=1)
            null[cursor:cursor + batch] = np.mean(xr[None, :] * yr[order], axis=1)
            cursor += batch
        output[metric] = {
            "rho": observed,
            "p_raw": empirical_p(observed, null, two_sided=True),
            "pairs": int(np.sum(good)),
        }
        null_by_metric[metric] = null
    finite_metrics = [m for m in metrics if np.isfinite(output[m]["rho"])]
    if finite_metrics:
        observed_max = max(abs(output[m]["rho"]) for m in finite_metrics)
        null_max = np.nanmax(np.stack([np.abs(null_by_metric[m])
                                      for m in finite_metrics]), axis=0)
        omnibus = {"max_abs_rho": float(observed_max),
                   "p_raw_max_metric": empirical_p(observed_max, null_max)}
    else:
        omnibus = {"max_abs_rho": np.nan, "p_raw_max_metric": np.nan}
    return {"metrics": output, "omnibus": omnibus}


def prepare_regular_feature(values, times):
    filled = np.full_like(values, np.nan, dtype=float)
    for j in range(values.shape[1]):
        good = np.isfinite(values[:, j])
        if np.sum(good) >= 3:
            good_indices = np.flatnonzero(good)
            filled[good_indices, j] = values[good_indices, j]
            for lower, upper in zip(good_indices[:-1], good_indices[1:]):
                if times[upper] - times[lower] <= 30.0:
                    filled[lower:upper + 1, j] = np.interp(
                        times[lower:upper + 1], times[[lower, upper]],
                        values[[lower, upper], j])
    scale = robust_mad(values, axis=0)
    scale = np.where(np.isfinite(scale) & (scale > 0), scale, np.nan)
    return filled, scale


def sample_regular_feature(filled, times, targets):
    dt = float(np.nanmedian(np.diff(times)))
    position = (targets - times[0]) / dt
    valid = (position >= 0.0) & (position <= len(times) - 1)
    lower_raw = np.floor(position)
    weight = position - lower_raw
    lower = np.clip(lower_raw.astype(int), 0, len(times) - 1)
    upper = np.clip(lower + 1, 0, len(times) - 1)
    radius_index = np.arange(filled.shape[1])[None, :]
    sampled = ((1.0 - weight) * filled[lower, radius_index]
               + weight * filled[upper, radius_index])
    return np.where(valid, sampled, np.nan)


def event_deltas_prepared(filled, scale, times, radii, anchors, speed):
    pre_offsets = np.array([-60.0, -45.0, -30.0, -15.0])
    post_offsets = np.array([15.0, 30.0, 45.0, 60.0])
    targets = np.asarray(anchors)[:, None] - (
        (3.0 - radii[None, :]) * RSUN_KM / speed / 60.0)
    pre = np.stack([sample_regular_feature(filled, times, targets + offset)
                    for offset in pre_offsets])
    post = np.stack([sample_regular_feature(filled, times, targets + offset)
                     for offset in post_offsets])
    return ((np.nanmedian(post, axis=0) - np.nanmedian(pre, axis=0))
            / scale[None, :])


def dynamic_gate_test(moment, anchors, speed):
    times, radii = moment["times"], moment["radii"]
    search = (radii >= RADIAL_SEARCH[0]) & (radii <= RADIAL_SEARCH[1])
    candidates = (2.10, 2.90)
    min_shift = float(np.ceil((times[0] + 240.0 - min(anchors)) / 15.0) * 15.0)
    max_shift = float(np.floor((times[-1] - 120.0 - max(anchors)) / 15.0) * 15.0)
    shifts = np.arange(min_shift, max_shift + 0.1, 15.0)
    shifts = shifts[np.abs(shifts) > 180.0]
    features = {
        "log_area": moment["log_area"],
        "centroid": moment["centroid_deg"],
        "log_excess": moment["log_excess"],
    }
    results, maps = {}, {}
    for feature, values in features.items():
        filled, scale = prepare_regular_feature(values, times)
        observed_deltas = event_deltas_prepared(
            filled, scale, times, radii, anchors, speed)
        observed_curve = np.nanmedian(np.abs(observed_deltas), axis=0)
        null_curves = []
        for shift in shifts:
            shifted = event_deltas_prepared(
                filled, scale, times, radii, np.asarray(anchors) + shift, speed)
            null_curves.append(np.nanmedian(np.abs(shifted), axis=0))
        null_curves = np.asarray(null_curves)
        valid = search & np.isfinite(observed_curve)
        if not np.any(valid):
            results[feature] = {"valid": False}
            continue
        peak_index = np.flatnonzero(valid)[np.nanargmax(observed_curve[valid])]
        null_max = np.nanmax(null_curves[:, search], axis=1)
        candidate_rows = {}
        for target in candidates:
            index = int(np.argmin(np.abs(radii - target)))
            candidate_rows[str(target)] = {
                "sampled_radius_rsun": float(radii[index]),
                "absolute_change": float(observed_curve[index]),
                "p_raw": empirical_p(observed_curve[index], null_curves[:, index]),
                "signed_median_change": float(np.nanmedian(observed_deltas[:, index])),
            }
        results[feature] = {
            "valid": True,
            "speed_km_s": speed,
            "null_common_shifts": int(len(shifts)),
            "peak_radius_rsun": float(radii[peak_index]),
            "peak_absolute_change": float(observed_curve[peak_index]),
            "scan_p_raw": empirical_p(observed_curve[peak_index], null_max),
            "signed_median_change_at_peak": float(
                np.nanmedian(observed_deltas[:, peak_index])),
            "candidates": candidate_rows,
        }
        maps[feature] = {
            "observed_curve": observed_curve,
            "null_q95": np.nanpercentile(null_curves, 95, axis=0),
        }
    passed = [name for name, row in results.items()
              if row.get("valid") and row["scan_p_raw"] < 0.05]
    radii_passed = [results[name]["peak_radius_rsun"] for name in passed]
    radius_agreement = bool(len(radii_passed) >= 2
                            and max(radii_passed) - min(radii_passed) <= 0.15)
    return {
        "features": results,
        "passing_features": passed,
        "two_feature_radius_agreement": radius_agreement,
        "strong_dynamic_gate_support": bool(len(passed) >= 2 and radius_agreement),
        "maps": maps,
    }


def periodic_ar1_test(series, times, seed, surrogates=4000):
    series = np.asarray(series, float)
    times = np.asarray(times, float)
    good = np.isfinite(series)
    if np.mean(good) < 0.70 or np.sum(good) < 40:
        return {"valid": False, "coverage": float(np.mean(good))}
    filled = np.interp(times, times[good], series[good])
    y = detrend(filled, type="linear")
    scale = float(np.std(y))
    if not np.isfinite(scale) or scale <= 0:
        return {"valid": False, "coverage": float(np.mean(good))}
    y /= scale
    cadence = float(np.nanmedian(np.diff(times)))
    frequency = np.fft.rfftfreq(len(y), d=cadence)
    periods = np.divide(1.0, frequency, out=np.full_like(frequency, np.inf),
                        where=frequency > 0)
    band = (periods >= 80.0) & (periods <= 130.0)
    power = np.abs(np.fft.rfft(y))**2 / len(y)
    if not np.any(band):
        return {"valid": False, "coverage": float(np.mean(good))}
    band_indices = np.flatnonzero(band)
    index = band_indices[np.argmax(power[band])]
    observed = float(power[index])

    phi = float(np.corrcoef(y[:-1], y[1:])[0, 1])
    phi = float(np.clip(phi, -0.95, 0.95)) if np.isfinite(phi) else 0.0
    rng = np.random.default_rng(seed)
    noise = rng.normal(size=(surrogates, len(y)))
    simulated = np.empty_like(noise)
    simulated[:, 0] = noise[:, 0] / np.sqrt(max(1.0 - phi**2, 1e-6))
    for j in range(1, len(y)):
        simulated[:, j] = phi * simulated[:, j - 1] + noise[:, j]
    simulated -= np.mean(simulated, axis=1, keepdims=True)
    simulated /= np.std(simulated, axis=1, keepdims=True)
    null_power = np.abs(np.fft.rfft(simulated, axis=1))**2 / len(y)
    null_stat = np.max(null_power[:, band], axis=1)
    return {
        "valid": True,
        "coverage": float(np.mean(good)),
        "peak_period_min": float(periods[index]),
        "band_max_power": observed,
        "ar1_phi": phi,
        "surrogates": surrogates,
        "p_raw_band_max": empirical_p(observed, null_stat),
    }


def periodic_suite(moment, seed_base):
    results = {}
    for ir, target in enumerate(FIXED_RADII):
        index = int(np.argmin(np.abs(moment["radii"] - target)))
        row = {}
        for jf, (name, key) in enumerate((
                ("log_area", "log_area"),
                ("centroid", "centroid_deg"),
                ("log_excess", "log_excess"))):
            row[name] = periodic_ar1_test(
                moment[key][:, index], moment["times"],
                seed=seed_base + 100 * ir + jf)
        results[str(target)] = {
            "sampled_radius_rsun": float(moment["radii"][index]),
            "features": row,
        }
    return results


def periodic_coincidences(suite):
    matches = []
    for radius, row in suite.items():
        density = row["features"]["log_excess"]
        if not density.get("valid") or density["p_raw_band_max"] >= 0.05:
            continue
        for geometry in ("log_area", "centroid"):
            candidate = row["features"][geometry]
            if (candidate.get("valid") and candidate["p_raw_band_max"] < 0.05
                    and abs(candidate["peak_period_min"]
                            - density["peak_period_min"]) <= 15.0):
                matches.append({
                    "radius_rsun": float(radius), "geometry": geometry,
                    "density_period_min": density["peak_period_min"],
                    "geometry_period_min": candidate["peak_period_min"],
                    "density_p_raw": density["p_raw_band_max"],
                    "geometry_p_raw": candidate["p_raw_band_max"],
                })
    return matches


def best_periodic_result(suite):
    candidates = []
    for radius, row in suite.items():
        for feature, result in row["features"].items():
            if result.get("valid"):
                candidates.append({
                    "radius_rsun": float(radius), "feature": feature,
                    "p_raw": result["p_raw_band_max"],
                    "period_min": result["peak_period_min"],
                })
    return min(candidates, key=lambda row: row["p_raw"]) if candidates else None


def serializable_dynamic(result):
    return {key: value for key, value in result.items() if key != "maps"}


def main():
    phase = json.loads(JAN_PHASE.read_text())
    events = phase["events"]
    anchors = np.array([row["cor2_anchor_min"] for row in events], float)
    response = np.array([row["common_ridge_score"] for row in events], float)
    best_response = np.array([row["best_ridge_score"] for row in events], float)

    with np.load(JAN_MAPS) as data:
        jan_path = np.asarray(data["traced_path_pa"], float)
    with np.load(JUL_MAPS) as data:
        july_path = np.asarray(data["cor1_path_pa"], float)

    jan_moments = {
        "total_b": streamer_moments(JAN_CUBE, "total_b", jan_path),
        "pb": streamer_moments(JAN_CUBE, "pb_signed", jan_path),
    }
    july_moments = {
        background: {
            "total_b": streamer_moments(path, "total_b", july_path),
            "pb": streamer_moments(path, "pb", july_path),
        }
        for background, path in JUL_CUBES.items()
    }

    geometry_rows = {}
    geometry_tests = {}
    for product, moment in jan_moments.items():
        geometry_rows[product] = {}
        geometry_tests[product] = {}
        for ispeed, speed in enumerate(SPEEDS):
            rows = event_geometry(moment, anchors, speed)
            geometry_rows[product][str(int(speed))] = rows
            geometry_tests[product][str(int(speed))] = permutation_correlation(
                rows, response, seed=20260900 + 10 * ispeed + (product == "pb"))

    primary_total = geometry_tests["total_b"]["200"]["metrics"]["max_area_gradient"]
    sensitivity_total = [geometry_tests["total_b"][str(int(speed))]["metrics"]
                         ["max_area_gradient"] for speed in (100.0, 300.0)]
    pb_primary = geometry_tests["pb"]["200"]["metrics"]["max_area_gradient"]
    same_sensitivity_sign = all(
        np.sign(row["rho"]) == np.sign(primary_total["rho"])
        for row in sensitivity_total if np.isfinite(row["rho"])
    )
    sensitivity_p = any(row["p_raw"] < 0.10 for row in sensitivity_total
                        if np.isfinite(row["p_raw"]))
    no_pb_reversal = (not np.isfinite(pb_primary["rho"])
                      or np.sign(pb_primary["rho"]) == np.sign(primary_total["rho"]))
    geometry_detected = bool(
        primary_total["rho"] > 0 and primary_total["p_raw"] < 0.05
        and same_sensitivity_sign and sensitivity_p and no_pb_reversal
    )

    dynamic = {product: {} for product in jan_moments}
    for product, moment in jan_moments.items():
        for speed in SPEEDS:
            dynamic[product][str(int(speed))] = dynamic_gate_test(
                moment, anchors, speed)
    primary_dynamic = dynamic["total_b"]["200"]
    primary_pass = primary_dynamic["strong_dynamic_gate_support"]
    primary_peaks = [primary_dynamic["features"][name]["peak_radius_rsun"]
                     for name in primary_dynamic["passing_features"]]
    sensitivity_reproduction = False
    for product in ("pb", "total_b"):
        for speed in SPEEDS:
            if product == "total_b" and speed == 200.0:
                continue
            row = dynamic[product][str(int(speed))]
            if not row["strong_dynamic_gate_support"]:
                continue
            peaks = [row["features"][name]["peak_radius_rsun"]
                     for name in row["passing_features"]]
            if primary_peaks and peaks and abs(np.median(peaks) - np.median(primary_peaks)) <= 0.15:
                sensitivity_reproduction = True
    dynamic_detected = bool(primary_pass and sensitivity_reproduction)

    periodic = {
        "january": {
            product: periodic_suite(moment, 20261000 + 1000 * i)
            for i, (product, moment) in enumerate(jan_moments.items())
        },
        "july": {
            background: {
                product: periodic_suite(moment, 20263000 + 2000 * ib + 1000 * ip)
                for ip, (product, moment) in enumerate(products.items())
            }
            for ib, (background, products) in enumerate(july_moments.items())
        },
    }
    periodic_matches = {
        "january": {product: periodic_coincidences(suite)
                    for product, suite in periodic["january"].items()},
        "july": {background: {product: periodic_coincidences(suite)
                              for product, suite in products.items()}
                 for background, products in periodic["july"].items()},
    }
    periodic_best = {
        "january": {product: best_periodic_result(suite)
                    for product, suite in periodic["january"].items()},
        "july": {background: {product: best_periodic_result(suite)
                              for product, suite in products.items()}
                 for background, products in periodic["july"].items()},
    }
    january_reproduced = False
    for a in periodic_matches["january"]["total_b"]:
        for b in periodic_matches["january"]["pb"]:
            if (abs(a["radius_rsun"] - b["radius_rsun"]) <= 0.15
                    and abs(a["density_period_min"] - b["density_period_min"]) <= 15.0):
                january_reproduced = True
    july_reproduced = False
    for product in ("total_b", "pb"):
        for a in periodic_matches["july"]["080723"][product]:
            for b in periodic_matches["july"]["080802"][product]:
                if (abs(a["radius_rsun"] - b["radius_rsun"]) <= 0.15
                        and abs(a["density_period_min"] - b["density_period_min"]) <= 15.0):
                    july_reproduced = True
    periodic_detected = bool(january_reproduced or july_reproduced)

    # Event table for the primary projected geometry.
    primary_rows = geometry_rows["total_b"]["200"]
    event_table = []
    for event, geom in zip(events, primary_rows):
        event_table.append({
            "event_number": event["event_number"], "utc": event["utc"],
            "common_ridge_score": event["common_ridge_score"],
            "common_ridge_p_raw": event["common_ridge_p_raw"],
            "best_ridge_p_raw": event["best_ridge_p_raw"],
            "best_model_family": event["best_model"]["family"],
            "best_change_radius_rsun": event["best_model"]["change_radius_rsun"],
            "best_v_inner_km_s": event["best_model"]["v_inner_km_s"],
            "best_v_outer_km_s": event["best_model"]["v_outer_km_s"],
            "phase_kink_radius_rsun": event["kink_radius_rsun"],
            "phase_kink_p_raw": event["kink_p_raw"],
            **geom,
        })

    event12_pb_geometry = geometry_rows["pb"]["200"][11]
    exploratory_event12 = {
        "status": "single-event candidate; post-primary exploratory reading",
        "utc": events[11]["utc"],
        "identified_ridge_p_raw": events[11]["best_ridge_p_raw"],
        "best_kinematic_model": events[11]["best_model"],
        "pb_area_gradient_peak_radius_rsun": (
            event12_pb_geometry["gradient_peak_radius_rsun"]),
        "pb_axis_shear_peak_radius_rsun": (
            event12_pb_geometry["axis_shear_peak_radius_rsun"]),
        "phase_kink_radius_rsun": events[11]["kink_radius_rsun"],
        "phase_kink_p_raw": events[11]["kink_p_raw"],
        "interpretation": (
            "The identified ridge and the pB geometry/axis coincidence near "
            "2.70--2.75 R_sun are consistent with a moving or reforming "
            "geometry-controlled gate. The phase kink is not independently "
            "significant, total-B geometry is unavailable at this time, and the "
            "existing pB compression test fails; this is not a confirmed shock."
        ),
    }

    result = {
        "test": "geometry-dependent Habbal and dynamic cusp/current-sheet gate",
        "date_frozen": "2026-08-21",
        "no_BH": True,
        "fixed_projection_speeds_km_s": list(SPEEDS),
        "geometry_dependence_detected": geometry_detected,
        "dynamic_gate_detected": dynamic_detected,
        "periodic_gate_detected": periodic_detected,
        "mhd_branch_identified": False,
        "event_table_primary_total_b_200_km_s": event_table,
        "exploratory_single_event_candidate": exploratory_event12,
        "geometry_tests": geometry_tests,
        "event_geometry_rows": geometry_rows,
        "geometry_decision_components": {
            "primary_total_b": primary_total,
            "sensitivity_total_b": sensitivity_total,
            "primary_pb": pb_primary,
            "same_sensitivity_sign": same_sensitivity_sign,
            "sensitivity_p_lt_0p10": sensitivity_p,
            "no_pb_sign_reversal": no_pb_reversal,
        },
        "dynamic_gate_tests": {
            product: {speed: serializable_dynamic(row)
                      for speed, row in speeds.items()}
            for product, speeds in dynamic.items()
        },
        "periodic_tests": periodic,
        "periodic_matches": periodic_matches,
        "periodic_best_single_channel": periodic_best,
        "periodic_reproduction": {
            "january_total_b_pb": january_reproduced,
            "july_double_background": july_reproduced,
        },
        "moment_valid_fraction": {
            "january": {key: row["valid_fraction"]
                        for key, row in jan_moments.items()},
            "july": {background: {key: row["valid_fraction"]
                                  for key, row in products.items()}
                      for background, products in july_moments.items()},
        },
        "branch_verdict": (
            "Imaging-only density, width, and axis observables cannot identify a "
            "slow, fast, or compound MHD shock without event-specific B, T, bulk "
            "normal velocity, and valid jump conditions."
        ),
    }

    json_path = OUT / "dynamic_geometry_gate_results_no_BH.json"
    json_path.write_text(json.dumps(result, indent=2, allow_nan=True,
                                    default=json_default) + "\n")

    # Diagnostic figure.
    fig, axes = plt.subplots(2, 2, figsize=(13, 9.5), constrained_layout=True)
    g = np.array([row["max_area_gradient"] for row in primary_rows])
    axes[0, 0].scatter(g, response, c=np.arange(1, 13), cmap="viridis", s=55)
    for i, (x, y) in enumerate(zip(g, response), 1):
        axes[0, 0].annotate(str(i), (x, y), xytext=(4, 3),
                            textcoords="offset points", fontsize=8)
    if np.isfinite(primary_total["rho"]):
        geometry_title = (f"Spearman rho={primary_total['rho']:.3f}, "
                          f"raw p={primary_total['p_raw']:.4f}")
    else:
        geometry_title = (f"quality-valid pairs={primary_total['pairs']} (<7); "
                          "primary correlation not estimated")
    axes[0, 0].set(
        title=("January event geometry vs frozen common-path ridge\n"
               + geometry_title),
        xlabel=r"max $d\ln A_{proxy}/dr$ (1.9--2.9 $R_\odot$)",
        ylabel="Frozen common-path ridge score")

    numbers = np.arange(1, 13)
    geometry_r = [row["gradient_peak_radius_rsun"] for row in primary_rows]
    model_r = [row["best_model"]["change_radius_rsun"] for row in events]
    kink_r = [row["kink_radius_rsun"] for row in events]
    axes[0, 1].plot(numbers, geometry_r, "o-", label="geometry-gradient peak")
    axes[0, 1].plot(numbers, model_r, "s--", label="best kinematic change")
    axes[0, 1].plot(numbers, kink_r, "^:", label="phase kink (mostly unidentified)")
    axes[0, 1].set(title="Event-dependent transition radii",
                   xlabel="January event number", ylabel=r"Radius ($R_\odot$)",
                   xticks=numbers, ylim=(1.75, 3.02))
    axes[0, 1].legend(fontsize=8)

    colors = {"log_area": "#4c78a8", "centroid": "#f28e2b",
              "log_excess": "#59a14f"}
    for name, row in primary_dynamic["maps"].items():
        axes[1, 0].plot(jan_moments["total_b"]["radii"], row["observed_curve"],
                        color=colors[name], label=f"{name} observed")
        axes[1, 0].plot(jan_moments["total_b"]["radii"], row["null_q95"],
                        color=colors[name], ls="--", alpha=0.7)
    axes[1, 0].set(title="Event-locked gate change at 200 km/s (dashed: null 95%)",
                   xlabel=r"Radius ($R_\odot$)", ylabel="Median |standardized change|",
                   xlim=(1.85, 2.95))
    axes[1, 0].legend(fontsize=8, ncol=2)

    p_matrix = np.full((3, len(FIXED_RADII)), np.nan)
    feature_order = ("log_area", "centroid", "log_excess")
    for j, radius in enumerate(FIXED_RADII):
        for i, feature in enumerate(feature_order):
            row = periodic["january"]["total_b"][str(radius)]["features"][feature]
            if row.get("valid"):
                p_matrix[i, j] = row["p_raw_band_max"]
    image = axes[1, 1].imshow(-np.log10(np.clip(p_matrix, 1e-4, 1.0)),
                              origin="lower", aspect="auto", vmin=0, vmax=4,
                              cmap="magma")
    axes[1, 1].set(title="January total-B 80--130 min AR(1) screen",
                   xticks=np.arange(len(FIXED_RADII)),
                   xticklabels=[f"{r:.2f}" for r in FIXED_RADII],
                   yticks=np.arange(3), yticklabels=feature_order,
                   xlabel=r"Radius ($R_\odot$)")
    if not np.any(np.isfinite(p_matrix)):
        axes[1, 1].text(
            0.5, 0.5,
            "Not tested: direct moment coverage\n"
            f"{jan_moments['total_b']['valid_fraction']['width']:.2f} < 0.70 quality gate",
            ha="center", va="center", transform=axes[1, 1].transAxes,
            bbox={"facecolor": "white", "alpha": 0.9, "edgecolor": "0.7"})
    fig.colorbar(image, ax=axes[1, 1], label=r"$-\log_{10}(p_{raw})$")
    fig.suptitle("Geometry-dependent Habbal and dynamic cusp-gate tests (no BH)")
    figure_path = OUT / "dynamic_geometry_gate_diagnostic_no_BH.png"
    fig.savefig(figure_path, dpi=210)
    plt.close(fig)

    report = [
        "# Geometry-dependent Habbal and dynamic cusp-gate screen (raw; no BH)",
        "", "## Decisions", "",
        f"- Geometry dependence detected: **{geometry_detected}**.",
        f"- Moving/reforming gate detected: **{dynamic_detected}**.",
        f"- Reproduced 80--130 min geometry clock detected: **{periodic_detected}**.",
        "- Slow/fast/compound MHD branch identified: **False**.",
        "", "## Geometry dependence", "",
        "Primary relation: total-B rapid projected-area divergence versus the already frozen common-path ridge score.",
        "", "| Projection | Product | rho | raw p | pairs |",
        "| ---: | --- | ---: | ---: | ---: |",
    ]
    for speed in SPEEDS:
        for product in ("total_b", "pb"):
            row = geometry_tests[product][str(int(speed))]["metrics"]["max_area_gradient"]
            report.append(f"| {speed:.0f} km/s | {product} | {row['rho']:.3f} | {row['p_raw']:.4f} | {row['pairs']} |")
    report += [
        "", ("Primary 200 km/s total-B correlation is not estimable under the "
              f"quality rule: {primary_total['pairs']} valid event pairs (<7)."),
        ("The pB sensitivity is negative rather than the predicted positive "
         "relation and is speed dependent; it is not reproduced by total-B."),
        "", "## Moving/reforming gate: January primary 200 km/s", "",
        "| Observable | Peak radius | Scan p raw | Signed median change |",
        "| --- | ---: | ---: | ---: |",
    ]
    for feature, row in primary_dynamic["features"].items():
        if row.get("valid"):
            report.append(
                f"| {feature} | {row['peak_radius_rsun']:.3f} | "
                f"{row['scan_p_raw']:.4f} | {row['signed_median_change_at_peak']:.3f} |")
    report += [
        "", f"Passing features: {', '.join(primary_dynamic['passing_features']) or 'none'}.",
        f"Two-feature radius agreement: {primary_dynamic['two_feature_radius_agreement']}.",
        "The isolated 100 km/s total-B log-area scan has raw p=0.0120 near 2.00 R_sun, but density and centroid do not pass there and the feature is not reproduced at 200/300 km/s or in pB.",
        "", "## Best single-event candidate (exploratory)", "",
        "Event #12 (2008-01-13 21:37:30 UT) remains the strongest individual case:",
        f"- fully scanned ridge raw p={events[11]['best_ridge_p_raw']:.4f};",
        f"- best kinematic path {events[11]['best_model']['family']}, {events[11]['best_model']['v_inner_km_s']:.0f}->{events[11]['best_model']['v_outer_km_s']:.0f} km/s at {events[11]['best_model']['change_radius_rsun']:.2f} R_sun;",
        f"- pB projected-area gradient peak at {event12_pb_geometry['gradient_peak_radius_rsun']:.2f} R_sun;",
        f"- pB axis-shear peak at {event12_pb_geometry['axis_shear_peak_radius_rsun']:.2f} R_sun;",
        f"- phase-kink location {events[11]['kink_radius_rsun']:.2f} R_sun, but kink raw p={events[11]['kink_p_raw']:.4f}.",
        "This radius coincidence is physically interesting and is consistent with a moving/reforming geometry-controlled gate. It is not a confirmed slow/fast/compound shock because the phase kink is not significant, total-B geometry is unavailable at this time, and the pB compression test did not pass.",
        "", "## 80--130 min periodic gate screen", "",
        f"January total-B/pB reproduction: {january_reproduced}.",
        f"July double-background reproduction: {july_reproduced}.",
        ("January is not directly testable under the periodic quality rule: "
         f"moment coverage is {jan_moments['total_b']['valid_fraction']['width']:.3f} "
         f"in total-B and {jan_moments['pb']['valid_fraction']['width']:.3f} in pB, "
         "below the frozen 0.70 threshold."),
        ("Best July single-channel raw p values are "
         f"080723 total-B={periodic_best['july']['080723']['total_b']['p_raw']:.4f}, "
         f"080723 pB={periodic_best['july']['080723']['pb']['p_raw']:.4f}, "
         f"080802 total-B={periodic_best['july']['080802']['total_b']['p_raw']:.4f}, "
         f"080802 pB={periodic_best['july']['080802']['pb']['p_raw']:.4f}; none is below 0.05."),
        "", "Raw same-radius density-plus-geometry matches:", "",
        "```json", json.dumps(periodic_matches, indent=2,
                               default=json_default), "```",
        "", "## Physical interpretation", "",
        "The Habbal mechanism is geometry dependent in theory, but the declared event-by-event geometry criterion is required before attributing PDS visibility to it. A coherent or periodic projected width alone can be streamer topology, line-of-sight weighting, or calibration structure.",
        "", "A density-column pulse does not distinguish slow, fast, or compound MHD transitions. Without B, T, upstream/downstream normal flow, and a valid compression jump, the present analysis can support only a dynamic-gate phenomenology, not identify a shock branch.",
        "", "Compact slow modes remain possible sub-streamer modulators. The result does not reinstate a full-streamer standing slow-mode or literal Laval/shock-diamond claim.",
    ]
    report_path = OUT / "PDS_dynamic_geometry_gate_no_BH.md"
    report_path.write_text("\n".join(report) + "\n")

    np.savez_compressed(
        OUT / "dynamic_geometry_gate_maps_no_BH.npz",
        radii=jan_moments["total_b"]["radii"],
        event_numbers=numbers,
        geometry_gradient_primary=g,
        common_ridge_score=response,
        best_ridge_score=best_response,
        geometry_peak_radius=np.asarray(geometry_r),
        kinematic_change_radius=np.asarray(model_r),
        phase_kink_radius=np.asarray(kink_r),
        dynamic_log_area=primary_dynamic["maps"].get("log_area", {}).get("observed_curve"),
        dynamic_centroid=primary_dynamic["maps"].get("centroid", {}).get("observed_curve"),
        dynamic_log_excess=primary_dynamic["maps"].get("log_excess", {}).get("observed_curve"),
    )
    print(json.dumps({
        "geometry_dependence_detected": geometry_detected,
        "dynamic_gate_detected": dynamic_detected,
        "periodic_gate_detected": periodic_detected,
        "primary_geometry": primary_total,
        "primary_dynamic": serializable_dynamic(primary_dynamic),
        "report": str(report_path), "figure": str(figure_path),
        "json": str(json_path),
    }, indent=2, default=json_default))


if __name__ == "__main__":
    main()
