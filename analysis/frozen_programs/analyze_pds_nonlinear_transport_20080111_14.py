#!/usr/bin/env python3
"""Nonlinear travel-time test for the 11--14 January 2008 PDS event.

The 12 strong COR2 ridges are fixed from the previous blind COR2 selection.
STEREO-A/COR1 total-brightness JP2 frames are downloaded from Helioviewer and
sampled in a polar sector around the independently fixed COR2 PA.  A common
event-centred ridge is then tested with three transport families:

1. constant pattern speed;
2. piecewise acceleration across a common radius;
3. piecewise deceleration across a common radius.

Model selection is performed on alternating training events and scored on the
held-out events.  The full selection is repeated after common time shifts of
the complete COR2 event train.  Reported p values are raw empirical p values;
no BH adjustment is applied.

The JP2 data are suitable for timing and morphology, not absolute pB
photometry.  A shock claim would require Level-1/pB and an independent Mach
number / Rankine--Hugoniot consistency test.
"""

from __future__ import annotations

import io
import json
import math
import re
import time
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
from scipy.ndimage import gaussian_filter1d, map_coordinates, uniform_filter1d


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results" / "pds_20080111_14_nonlinear"
OUT.mkdir(parents=True, exist_ok=True)
CACHE = OUT / "cor1a_hv_polar_sector_15min.npz"

SOURCE_ID_COR1_A = 28
RSUN_KM = 695700.0
T0 = datetime(2008, 1, 11, tzinfo=timezone.utc)
START = datetime(2008, 1, 10, 12, 5, tzinfo=timezone.utc)
STOP = datetime(2008, 1, 14, 5, 50, tzinfo=timezone.utc)
CADENCE_MIN = 15.0
SEED_PA_DEG = 174.5

# Independently selected in COR2 at 3 R_sun.  These are the z >= 1.25 subset
# from pds_ridge_onset_results.json; they are not re-selected in COR1.
STRONG_EVENTS = (
    ("2008-01-11T14:37:30Z", 1.9045678101675219),
    ("2008-01-12T05:07:30Z", 1.4185207800009825),
    ("2008-01-12T11:37:30Z", 1.3572946759855757),
    ("2008-01-12T22:07:30Z", 2.0268880649347003),
    ("2008-01-13T01:07:30Z", 2.0289245855915374),
    ("2008-01-13T05:37:30Z", 1.2689771940366317),
    ("2008-01-13T10:07:30Z", 1.5508401029178402),
    ("2008-01-13T14:37:30Z", 1.5048374337932233),
    ("2008-01-13T15:37:30Z", 1.6497002418339775),
    ("2008-01-13T17:07:30Z", 1.6386389862359545),
    ("2008-01-13T19:37:30Z", 1.4046445783936397),
    ("2008-01-13T21:37:30Z", 1.9016556531675781),
)


def parse_utc(text: str) -> datetime:
    return datetime.fromisoformat(text.replace("Z", "+00:00"))


EVENT_TIMES = np.array(
    [(parse_utc(t) - T0).total_seconds() / 60.0 for t, _ in STRONG_EVENTS],
    dtype=float,
)
EVENT_STRENGTHS = np.array([z for _, z in STRONG_EVENTS], dtype=float)


def xml_value(blob: bytes, tag: str, default: float | str | None = None):
    pattern = rb"<" + re.escape(tag.encode()) + rb">(.*?)</" + re.escape(tag.encode()) + rb">"
    match = re.search(pattern, blob[:40000], flags=re.DOTALL)
    if not match:
        return default
    return match.group(1).decode("utf-8", "replace").strip()


def download_one(requested: datetime, retries: int = 4):
    query = urllib.parse.urlencode(
        {"date": requested.strftime("%Y-%m-%dT%H:%M:%SZ"), "sourceId": SOURCE_ID_COR1_A}
    )
    url = f"https://api.helioviewer.org/v2/getJP2Image/?{query}"
    error = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "PDS-nonlinear-test/1.0"})
            with urllib.request.urlopen(req, timeout=60) as response:
                blob = response.read()
            date_text = xml_value(blob, "DATE_OBS")
            if not date_text:
                raise RuntimeError("DATE_OBS absent from JP2 metadata")
            observed = datetime.fromisoformat(str(date_text)).replace(tzinfo=timezone.utc)
            image = np.asarray(Image.open(io.BytesIO(blob)).convert("L"), dtype=np.float32)
            crpix1 = float(xml_value(blob, "CRPIX1")) - 1.0
            crpix2 = float(xml_value(blob, "CRPIX2")) - 1.0
            rsun_arcsec = float(xml_value(blob, "RSUN"))
            cdelt = abs(float(xml_value(blob, "CDELT1")))
            return requested, observed, image, crpix1, crpix2, rsun_arcsec / cdelt
        except Exception as exc:  # network retries are intentionally narrow
            error = exc
            time.sleep(1.5 * (attempt + 1))
    raise RuntimeError(f"Failed {requested.isoformat()}: {error}")


def polar_sample(image, cx, cy, rsun_px, radii, pa_deg):
    theta = np.deg2rad(pa_deg)[None, :]
    rr = radii[:, None] * rsun_px
    xx = cx + rr * np.sin(theta)
    yy = cy - rr * np.cos(theta)
    return map_coordinates(image, [yy, xx], order=1, mode="constant", cval=np.nan)


def requested_times():
    out = []
    current = START
    while current <= STOP:
        out.append(current)
        current += timedelta(minutes=CADENCE_MIN)
    return out


def build_or_load_cube():
    if CACHE.exists():
        data = np.load(CACHE)
        return {key: data[key] for key in data.files}

    radii = np.arange(1.45, 3.351, 0.025, dtype=float)
    pa = np.arange(140.0, 210.001, 1.0, dtype=float)
    requests = requested_times()
    results = []
    print(f"Downloading {len(requests)} COR1-A JP2 frames...", flush=True)
    with ThreadPoolExecutor(max_workers=8) as pool:
        future_map = {pool.submit(download_one, dt): dt for dt in requests}
        for number, future in enumerate(as_completed(future_map), 1):
            results.append(future.result())
            if number % 24 == 0 or number == len(requests):
                print(f"  COR1 frames: {number}/{len(requests)}", flush=True)

    # Nearest-image requests can duplicate a frame across a data gap.  Keep one
    # copy per DATE_OBS, reject a match farther than half a requested cadence.
    unique = {}
    for requested, observed, image, cx, cy, rsun_px in results:
        delta = abs((observed - requested).total_seconds()) / 60.0
        if delta <= CADENCE_MIN / 2 + 0.6:
            unique[observed] = (image, cx, cy, rsun_px)

    times = sorted(unique)
    cube = np.empty((len(times), len(radii), len(pa)), dtype=np.float32)
    for i, observed in enumerate(times):
        image, cx, cy, rsun_px = unique[observed]
        cube[i] = polar_sample(image, cx, cy, rsun_px, radii, pa)
    minutes = np.array([(dt - T0).total_seconds() / 60.0 for dt in times], dtype=float)
    np.savez_compressed(
        CACHE,
        intensity=cube,
        minutes=minutes,
        radii=radii,
        pa=pa,
        dates=np.array([dt.isoformat() for dt in times]),
        source=np.array("Helioviewer JP2 sourceId=28; 8-bit processed COR1-A total brightness"),
    )
    return {"intensity": cube, "minutes": minutes, "radii": radii, "pa": pa}


def regularize(minutes, cube, cadence=CADENCE_MIN):
    grid = np.arange(
        math.ceil(np.nanmin(minutes) / cadence) * cadence,
        math.floor(np.nanmax(minutes) / cadence) * cadence + 0.1,
        cadence,
    )
    flat = np.asarray(cube, float).reshape(len(minutes), -1)
    out = np.empty((len(grid), flat.shape[1]), dtype=np.float64)
    for j in range(flat.shape[1]):
        good = np.isfinite(flat[:, j])
        if np.sum(good) >= 2:
            out[:, j] = np.interp(grid, minutes[good], flat[good, j])
        else:
            out[:, j] = np.nan
    return grid, out.reshape((len(grid),) + cube.shape[1:])


def trace_path(cube, radii, pa, seed_pa=SEED_PA_DEG):
    """Trace the static southern ray while limiting radial PA wandering."""
    static = np.nanmedian(cube, axis=0)
    static = gaussian_filter1d(static, 1.4, axis=1, mode="nearest")
    path = np.full(len(radii), seed_pa, dtype=float)
    anchor = int(np.argmin(np.abs(radii - 3.0)))

    def pick(j, center, half_window):
        use = np.abs(pa - center) <= half_window
        profile = static[j].copy()
        profile[~use] = -np.inf
        return float(pa[int(np.nanargmax(profile))]) if np.any(np.isfinite(profile)) else center

    path[anchor] = pick(anchor, seed_pa, 10.0)
    for j in range(anchor - 1, -1, -1):
        path[j] = pick(j, path[j + 1], 4.0)
    for j in range(anchor + 1, len(radii)):
        path[j] = pick(j, path[j - 1], 4.0)
    path = gaussian_filter1d(path, 1.6, mode="nearest")
    return path, static


def sample_path(cube, pa, path, half_width=2.5):
    out = np.empty((cube.shape[0], cube.shape[1]), dtype=float)
    for j, center in enumerate(path):
        use = np.abs(pa - center) <= half_width
        out[:, j] = np.nanmean(cube[:, j, use], axis=1)
    return out


def bff_maps(values, cadence=CADENCE_MIN, short_minutes=75.0, long_minutes=618.0):
    x = np.asarray(values, dtype=float)
    s_short = max(short_minutes / cadence / 2.355, 0.75)
    s_long = max(long_minutes / cadence / 2.355, 1.5)
    raw = gaussian_filter1d(x, s_short, axis=0, mode="nearest") - gaussian_filter1d(
        x, s_long, axis=0, mode="nearest"
    )
    center = np.nanmedian(raw, axis=0)
    mad = 1.4826 * np.nanmedian(np.abs(raw - center[None]), axis=0)
    std = np.nanstd(raw, axis=0)
    scale = np.where(mad > 1e-8, mad, np.where(std > 1e-8, std, 1.0))
    z = np.clip((raw - center[None]) / scale[None], -6.0, 6.0)
    z = uniform_filter1d(z, 13, axis=1, mode="nearest")
    return raw, z


def event_stack(grid, zmap, anchors, lags):
    stack = np.empty((len(lags), zmap.shape[1]), dtype=float)
    targets = anchors[:, None] + lags[None, :]
    for j in range(zmap.shape[1]):
        samples = np.array(
            [np.interp(targets[k], grid, zmap[:, j], left=np.nan, right=np.nan)
             for k in range(len(anchors))]
        )
        stack[:, j] = np.nanmean(samples, axis=0)
    return stack


@dataclass(frozen=True)
class Model:
    family: str
    v_inner: float
    v_outer: float
    change_radius: float | None
    tau: np.ndarray


def travel_time_piecewise(radii, v_inner, v_outer, change_radius):
    radii = np.asarray(radii, float)
    tau = np.empty_like(radii)
    outer = radii >= change_radius
    tau[outer] = -(3.0 - radii[outer]) * RSUN_KM / v_outer / 60.0
    tau[~outer] = -(
        (3.0 - change_radius) * RSUN_KM / v_outer
        + (change_radius - radii[~outer]) * RSUN_KM / v_inner
    ) / 60.0
    return tau


def model_families(radii):
    speeds = np.array([20, 25, 30, 40, 50, 60, 80, 100, 125, 150, 175, 200, 250, 300], float)
    constants = [
        Model("constant", v, v, None, -(3.0 - radii) * RSUN_KM / v / 60.0)
        for v in speeds
    ]
    changes = np.arange(1.8, 2.901, 0.1)
    acceleration, deceleration = [], []
    for rc in changes:
        for vin in speeds:
            for vout in speeds:
                if vout >= 1.25 * vin:
                    tau = travel_time_piecewise(radii, vin, vout, rc)
                    if np.nanmin(tau) >= -920:
                        acceleration.append(Model("acceleration", vin, vout, float(rc), tau))
                if vout <= 0.80 * vin:
                    tau = travel_time_piecewise(radii, vin, vout, rc)
                    if np.nanmin(tau) >= -920:
                        deceleration.append(Model("deceleration", vin, vout, float(rc), tau))
    return {"constant": constants, "acceleration": acceleration, "deceleration": deceleration}


def sample_models(stack, lags, models):
    """Return mean stacked ridge score for every precomputed model."""
    dt = float(np.median(np.diff(lags)))
    tau = np.stack([m.tau for m in models])
    pos = (tau - lags[0]) / dt
    lo = np.floor(pos).astype(int)
    hi = np.clip(lo + 1, 0, len(lags) - 1)
    lo = np.clip(lo, 0, len(lags) - 1)
    weight = pos - np.floor(pos)
    ridx = np.arange(stack.shape[1])[None, :]
    values = (1.0 - weight) * stack[lo, ridx] + weight * stack[hi, ridx]
    return np.nanmean(values, axis=1)


def score_fixed_model(stack, lags, model):
    return float(sample_models(stack, lags, [model])[0])


def cross_validated_at_shift(grid, zmap, shifted_events, lags, families):
    even = np.arange(len(shifted_events)) % 2 == 0
    folds = ((even, ~even), (~even, even))
    output = {}
    for name, models in families.items():
        held_out, chosen = [], []
        for train, valid in folds:
            train_stack = event_stack(grid, zmap, shifted_events[train], lags)
            train_scores = sample_models(train_stack, lags, models)
            best_index = int(np.nanargmax(train_scores))
            valid_stack = event_stack(grid, zmap, shifted_events[valid], lags)
            held_out.append(score_fixed_model(valid_stack, lags, models[best_index]))
            chosen.append(best_index)
        output[name] = {"cv_score": float(np.mean(held_out)), "fold_scores": held_out,
                        "chosen_indices": chosen}
    return output


def model_to_dict(model):
    return {
        "family": model.family,
        "v_inner_km_s": model.v_inner,
        "v_outer_km_s": model.v_outer,
        "change_radius_rsun": model.change_radius,
        "travel_time_at_1p55_min": float(np.interp(1.55, RADII_WORK, model.tau)),
    }


def empirical_tests(grid, zmap, lags, families):
    observed = cross_validated_at_shift(grid, zmap, EVENT_TIMES, lags, families)
    # A valid common shift must keep the slowest allowed precursor and all
    # anchor times inside the measured time range.
    min_lag = min(float(np.min(m.tau)) for group in families.values() for m in group)
    shift_lo = grid[0] - np.min(EVENT_TIMES) - min_lag
    shift_hi = grid[-1] - np.max(EVENT_TIMES)
    shifts = np.arange(math.ceil(shift_lo / CADENCE_MIN) * CADENCE_MIN,
                       math.floor(shift_hi / CADENCE_MIN) * CADENCE_MIN + 0.1,
                       CADENCE_MIN)
    null_shifts = shifts[np.abs(shifts) > 1e-6]
    null = {name: [] for name in families}
    null_improvements = {"acceleration_minus_constant": [], "deceleration_minus_constant": []}
    for number, shift in enumerate(null_shifts, 1):
        result = cross_validated_at_shift(grid, zmap, EVENT_TIMES + shift, lags, families)
        for name in families:
            null[name].append(result[name]["cv_score"])
        null_improvements["acceleration_minus_constant"].append(
            result["acceleration"]["cv_score"] - result["constant"]["cv_score"]
        )
        null_improvements["deceleration_minus_constant"].append(
            result["deceleration"]["cv_score"] - result["constant"]["cv_score"]
        )
        if number % 30 == 0:
            print(f"  null shifts: {number}/{len(null_shifts)}", flush=True)

    for name in families:
        values = np.asarray(null[name], float)
        obs = observed[name]["cv_score"]
        observed[name]["raw_p_common_shift"] = float((np.sum(values >= obs) + 1) / (len(values) + 1))
        observed[name]["null_mean"] = float(np.mean(values))
        observed[name]["null_std"] = float(np.std(values))
        observed[name]["null_scores"] = values.tolist()
        observed[name]["chosen_models"] = [
            model_to_dict(families[name][i]) for i in observed[name]["chosen_indices"]
        ]

    for key, values in null_improvements.items():
        family = key.split("_minus_")[0]
        obs = observed[family]["cv_score"] - observed["constant"]["cv_score"]
        values = np.asarray(values, float)
        observed[key] = {
            "observed_score_gain": float(obs),
            "raw_p_common_shift": float((np.sum(values >= obs) + 1) / (len(values) + 1)),
            "null_mean": float(np.mean(values)),
        }
    observed["null_shift_count"] = int(len(null_shifts))
    observed["null_shift_range_min"] = [float(np.min(null_shifts)), float(np.max(null_shifts))]
    return observed


def best_all_data_models(stack, lags, families):
    chosen = {}
    for name, models in families.items():
        scores = sample_models(stack, lags, models)
        index = int(np.nanargmax(scores))
        chosen[name] = (models[index], float(scores[index]))
    return chosen


def profile_diagnostics(grid, zmap, radii, model):
    amplitudes = []
    jitter = []
    event_offsets = []
    window = np.arange(-45.0, 45.1, CADENCE_MIN)
    for j, (r, tau) in enumerate(zip(radii, model.tau)):
        vals = np.array([
            [np.interp(event + tau + delta, grid, zmap[:, j]) for delta in window]
            for event in EVENT_TIMES
        ])
        stack = np.mean(vals, axis=0)
        amplitudes.append(float(np.interp(0.0, window, stack)))
        offsets = window[np.argmax(vals, axis=1)]
        event_offsets.append(offsets)
        jitter.append(float(1.4826 * np.median(np.abs(offsets - np.median(offsets)))))
    return np.asarray(amplitudes), np.asarray(jitter), np.asarray(event_offsets).T


def make_figure(grid, cube, radii_all, pa, path, zmap, lags, stack, families, best, tests):
    fig, axes = plt.subplots(2, 3, figsize=(16, 10), constrained_layout=True)
    static = np.nanmedian(cube, axis=0)
    vmax = np.nanpercentile(static, 99)
    axes[0, 0].pcolormesh(pa, radii_all, static, shading="auto", cmap="gray", vmin=0, vmax=vmax)
    axes[0, 0].plot(path, radii_all, color="#e45756", lw=2)
    axes[0, 0].axhline(3.0, color="cyan", ls="--", lw=1)
    axes[0, 0].set(title="COR1-A static polar sector and traced ray", xlabel="PA (deg)",
                   ylabel=r"Radius ($R_\odot$)")

    lim = np.nanpercentile(np.abs(stack), 98)
    axes[0, 1].pcolormesh(lags, RADII_WORK, stack.T, shading="auto", cmap="RdBu_r",
                          vmin=-lim, vmax=lim)
    colors = {"constant": "white", "acceleration": "lime", "deceleration": "gold"}
    for name, (model, _) in best.items():
        axes[0, 1].plot(model.tau, RADII_WORK, color=colors[name], lw=2,
                        label=f"{name}: {model.v_inner:.0f}/{model.v_outer:.0f} km s$^{{-1}}$")
    axes[0, 1].set(title="Stack of 12 COR2-selected events", xlabel=r"Minutes relative to 3 R$_\odot$ anchor",
                   ylabel=r"Radius ($R_\odot$)")
    axes[0, 1].legend(fontsize=8, loc="lower right")

    names = ["constant", "acceleration", "deceleration"]
    pvals = [tests[name]["raw_p_common_shift"] for name in names]
    bars = axes[0, 2].bar(names, -np.log10(pvals), color=[colors[n] for n in names], edgecolor="k")
    axes[0, 2].axhline(-np.log10(0.05), color="k", ls="--")
    axes[0, 2].set(title="Held-out common-shift test (no BH)", ylabel=r"$-\log_{10} p_{raw}$")
    for bar, p in zip(bars, pvals):
        axes[0, 2].text(bar.get_x() + bar.get_width()/2, bar.get_height()+0.03, f"p={p:.3f}",
                        ha="center", va="bottom", fontsize=9)

    for name, (model, _) in best.items():
        speed = np.where(RADII_WORK < (model.change_radius or -99), model.v_inner, model.v_outer)
        axes[1, 0].step(RADII_WORK, speed, where="mid", color=colors[name], lw=2, label=name)
    axes[1, 0].set(title="Best all-event kinematic profiles", xlabel=r"Radius ($R_\odot$)",
                   ylabel="Pattern speed (km s$^{-1}$)")
    axes[1, 0].legend()

    # Profile diagnostics follow the best held-out family by p value, while
    # retaining the all-data parameter fit only for visualization.
    chosen_family = min(names, key=lambda n: tests[n]["raw_p_common_shift"])
    chosen_model = best[chosen_family][0]
    amplitude, jitter, event_offsets = profile_diagnostics(grid, zmap, RADII_WORK, chosen_model)
    axes[1, 1].plot(RADII_WORK, amplitude, color="#2a9d8f", lw=2)
    if chosen_model.change_radius is not None:
        axes[1, 1].axvline(chosen_model.change_radius, color=colors[chosen_family], ls="--")
    axes[1, 1].axhline(0, color="0.5", lw=1)
    axes[1, 1].set(title=f"Stacked ridge amplitude along {chosen_family} path",
                   xlabel=r"Radius ($R_\odot$)", ylabel="Robust BFF z")

    axes[1, 2].plot(RADII_WORK, jitter, color="#6a4c93", lw=2)
    if chosen_model.change_radius is not None:
        axes[1, 2].axvline(chosen_model.change_radius, color=colors[chosen_family], ls="--")
    axes[1, 2].set(title="Event-to-event timing jitter (diagnostic)",
                   xlabel=r"Radius ($R_\odot$)", ylabel="Robust timing scatter (min)")
    fig.suptitle("11–14 January 2008: nonlinear COR1 travel-time test", fontsize=16)
    fig.savefig(OUT / "pds_20080111_14_nonlinear_transport_diagnostic.png", dpi=180)
    plt.close(fig)
    return chosen_family, chosen_model, amplitude, jitter, event_offsets


def write_report(tests, best, chosen_family, chosen_model, amplitude, jitter, path, frame_count):
    def fmt_model(model):
        if model.change_radius is None:
            return f"v={model.v_inner:.0f} km/s"
        return (f"v_inner={model.v_inner:.0f} km/s, v_outer={model.v_outer:.0f} km/s, "
                f"r_change={model.change_radius:.2f} R_sun")

    accel_models = tests["acceleration"]["chosen_models"]
    decel_models = tests["deceleration"]["chosen_models"]
    accel_stable = abs(accel_models[0]["change_radius_rsun"] - accel_models[1]["change_radius_rsun"]) <= 0.2
    decel_stable = abs(decel_models[0]["change_radius_rsun"] - decel_models[1]["change_radius_rsun"]) <= 0.2
    lines = [
        "# Nonlinear PDS transport test: 11–14 January 2008",
        "",
        "## Question",
        "",
        "Can the earlier failure of a fixed 200 km/s continuation be explained by a common",
        "acceleration or shock-like change in propagation speed inside COR1?",
        "",
        "## Frozen inputs",
        "",
        "- Event: 11–14 January 2008.",
        "- 12 strong COR2 ridges (z >= 1.25), selected previously at 3 R_sun.",
        f"- Fixed COR2 PA: {SEED_PA_DEG:.1f} deg; COR1 ray is traced around that PA.",
        f"- COR1-A frames used: {frame_count}, 15 min sampling.",
        "- No event was added, deleted, or moved after looking at COR1.",
        "- No BH correction is applied; all p values below are raw empirical p values.",
        "",
        "## Data caveat",
        "",
        "The timing map uses processed 8-bit Helioviewer COR1-A JP2 total-brightness images.",
        "Their embedded metadata show that the COR1 product was formed from the three",
        "polarization exposures, but JP2 values are not absolute pB photometry.  Therefore",
        "this is a timing/morphology test; compression ratios and Mach numbers require Level-1 pB.",
        "",
        "## Test",
        "",
        "Alternating halves of the 12-event train are used for two-fold cross-validation.",
        "For each fold, the training events select a constant-speed, accelerating, or",
        "decelerating piecewise travel-time path.  That frozen path is scored on the other",
        "six events.  Every speed and change-radius scan is repeated after common shifts of",
        "the full COR2 train, preserving its cadence and clustering.",
        "",
        "## Raw results",
        "",
        "| Model family | Held-out score | Raw common-shift p | Fold-selected models |",
        "| --- | ---: | ---: | --- |",
    ]
    for name in ("constant", "acceleration", "deceleration"):
        selected = tests[name]["chosen_models"]
        model_text = "; ".join(
            (f"{m['v_inner_km_s']:.0f} km/s" if m["change_radius_rsun"] is None else
             f"{m['v_inner_km_s']:.0f}->{m['v_outer_km_s']:.0f} km/s at {m['change_radius_rsun']:.2f} R_sun")
            for m in selected
        )
        lines.append(
            f"| {name} | {tests[name]['cv_score']:.4f} | {tests[name]['raw_p_common_shift']:.4f} | {model_text} |"
        )
    lines += [
        "",
        "### Does a speed break improve the held-out prediction?",
        "",
        f"- Acceleration minus constant: score gain={tests['acceleration_minus_constant']['observed_score_gain']:.4f}, "
        f"raw p={tests['acceleration_minus_constant']['raw_p_common_shift']:.4f}.",
        f"- Deceleration minus constant: score gain={tests['deceleration_minus_constant']['observed_score_gain']:.4f}, "
        f"raw p={tests['deceleration_minus_constant']['raw_p_common_shift']:.4f}.",
        f"- Acceleration change radius stable between folds within 0.2 R_sun: {accel_stable}.",
        f"- Deceleration change radius stable between folds within 0.2 R_sun: {decel_stable}.",
        "",
        "## All-event visualization fits (not used for p values)",
        "",
    ]
    for name, (model, score) in best.items():
        lines.append(f"- {name}: {fmt_model(model)}; stacked score={score:.4f}.")
    lines += [
        "",
        "## Physical decision rule",
        "",
        "A stationary acceleration/shock candidate requires all of the following:",
        "",
        "1. raw common-shift p < 0.05 on held-out events;",
        "2. a significant held-out gain relative to constant speed;",
        "3. a similar change radius in both folds;",
        "4. neighbouring-height ridge continuity and a repeatable amplitude/shape change;",
        "5. Level-1 pB plus an independent magnetosonic Mach-number check before calling it a shock.",
        "",
        "## Interpretation",
        "",
    ]
    candidate = (
        tests[chosen_family]["raw_p_common_shift"] < 0.05
        and chosen_family != "constant"
        and tests[f"{chosen_family}_minus_constant"]["raw_p_common_shift"] < 0.05
        and ((chosen_family == "acceleration" and accel_stable)
             or (chosen_family == "deceleration" and decel_stable))
    )
    if candidate:
        lines += [
            f"The held-out timing test supports a repeatable {chosen_family} break.  The",
            "result is a kinematic transition candidate, not yet a shock detection.  The",
            "next confirmation must use Level-1 pB to test compression and M_f > 1.",
        ]
    else:
        lines += [
            "The strict held-out criteria do not establish a stationary acceleration or",
            "shock-like break.  A visually improved all-event path is insufficient if the",
            "gain over constant speed or the change radius is not reproducible between folds.",
            "The high-COR1 aggregate onset can therefore still represent local generation,",
            "a time-dependent/moving boundary, or a visibility threshold.",
        ]
    lines += [
        "",
        "The fixed 200 km/s failure remains non-decisive: this test asks the stronger and",
        "physically relevant question of whether one nonlinear travel-time law predicts",
        "independent ridges.",
    ]
    (OUT / "PDS_20080111_14_nonlinear_transport_no_BH.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )
    return candidate


def main():
    global RADII_WORK
    data = build_or_load_cube()
    minutes, cube = regularize(data["minutes"], data["intensity"], CADENCE_MIN)
    # Remove frame-wise radial/common display scaling before following the ray.
    contrast = cube - np.nanmedian(cube, axis=2, keepdims=True)
    path, static = trace_path(cube, data["radii"], data["pa"])
    path_values = sample_path(contrast, data["pa"], path, half_width=2.5)
    raw_bff, z_all = bff_maps(path_values)
    use = (data["radii"] >= 1.55) & (data["radii"] <= 3.001)
    RADII_WORK = data["radii"][use]
    zmap = z_all[:, use]
    raw_map = raw_bff[:, use]

    lags = np.arange(-920.0, 0.1, CADENCE_MIN)
    families = model_families(RADII_WORK)
    print("Model counts:", {k: len(v) for k, v in families.items()}, flush=True)
    tests = empirical_tests(minutes, zmap, lags, families)
    stack = event_stack(minutes, zmap, EVENT_TIMES, lags)
    best = best_all_data_models(stack, lags, families)
    chosen_family, chosen_model, amplitude, jitter, event_offsets = make_figure(
        minutes, cube, data["radii"], data["pa"], path, zmap, lags, stack,
        families, best, tests,
    )
    candidate = write_report(
        tests, best, chosen_family, chosen_model, amplitude, jitter, path, len(data["minutes"])
    )

    output = {
        "event": "2008-01-11--14",
        "cor2_fixed_pa_deg": SEED_PA_DEG,
        "strong_event_count": int(len(EVENT_TIMES)),
        "strong_event_times_utc": [t for t, _ in STRONG_EVENTS],
        "strong_event_z": EVENT_STRENGTHS.tolist(),
        "cor1_frame_count": int(len(data["minutes"])),
        "cor1_time_range_min_from_20080111": [float(minutes[0]), float(minutes[-1])],
        "traced_path_pa_range_deg": [float(np.min(path)), float(np.max(path))],
        "model_counts": {k: len(v) for k, v in families.items()},
        "tests": tests,
        "best_all_event_models": {
            name: {**model_to_dict(model), "stacked_score": score}
            for name, (model, score) in best.items()
        },
        "diagnostic_selected_family": chosen_family,
        "diagnostic_selected_model": model_to_dict(chosen_model),
        "ridge_amplitude_profile": amplitude.tolist(),
        "timing_jitter_profile_min": jitter.tolist(),
        "event_timing_offsets_min": event_offsets.tolist(),
        "radii_rsun": RADII_WORK.tolist(),
        "strict_transition_candidate": bool(candidate),
        "no_BH": True,
        "data_caveat": "Helioviewer 8-bit processed COR1-A JP2; timing/morphology only, not absolute pB",
    }
    (OUT / "pds_20080111_14_nonlinear_transport_results.json").write_text(
        json.dumps(output, indent=2), encoding="utf-8"
    )
    np.savez_compressed(
        OUT / "pds_20080111_14_nonlinear_transport_maps.npz",
        minutes=minutes,
        radii=RADII_WORK,
        zmap=zmap.astype(np.float32),
        raw_bff=raw_map.astype(np.float32),
        lags=lags,
        event_stack=stack.astype(np.float32),
        traced_path_pa=path,
        event_times=EVENT_TIMES,
    )
    print(json.dumps({
        "p": {name: tests[name]["raw_p_common_shift"]
              for name in ("constant", "acceleration", "deceleration")},
        "gain_p": {
            "acceleration": tests["acceleration_minus_constant"]["raw_p_common_shift"],
            "deceleration": tests["deceleration_minus_constant"]["raw_p_common_shift"],
        },
        "best": {name: model_to_dict(model) for name, (model, _) in best.items()},
        "strict_transition_candidate": candidate,
    }, indent=2))


if __name__ == "__main__":
    main()
