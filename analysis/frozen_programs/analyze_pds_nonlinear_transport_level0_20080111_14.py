#!/usr/bin/env python3
"""Level-0 COR1 confirmation of the nonlinear PDS travel-time test.

This script downloads the three STEREO-A/COR1 polarizer exposures for each
15-minute sequence, subtracts the header bias, divides by exposure time, and
fits total brightness plus the tangential polarized component directly in a
solar polar sector.  Only the small polar arrays are retained.

The 12 COR2 anchors, PA, model families and train/validation split are frozen
before COR1 is inspected.  Product choice (tB or pB), speed and change radius
are selected on training ridges and scored on held-out ridges.  The complete
selection is repeated in a common-shift empirical null.  No BH adjustment is
applied.
"""

from __future__ import annotations

import io
import json
import math
import re
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import timedelta
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from scipy.ndimage import gaussian_filter1d, map_coordinates

import analyze_pds_nonlinear_transport_20080111_14 as core


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results" / "pds_20080111_14_nonlinear_level0"
OUT.mkdir(parents=True, exist_ok=True)
CACHE = OUT / "cor1a_level0_fitpol_sector_15min.npz"
TRIPLET_CACHE = OUT / "triplet_cache"
TRIPLET_CACHE.mkdir(parents=True, exist_ok=True)


def fetch(url, retries=4):
    error = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "PDS-Level0-test/1.0"})
            with urllib.request.urlopen(req, timeout=90) as response:
                return response.read()
        except Exception as exc:
            error = exc
            time.sleep(1.5 * (attempt + 1))
    raise RuntimeError(f"Download failed: {url}: {error}")


def read_image(blob):
    with fits.open(io.BytesIO(blob), memmap=False) as hdul:
        data = hdul[0].data.astype(np.float64)
        header = hdul[0].header.copy()
    bias = float(header.get("BIASMEAN", 0.0))
    exptime = float(header.get("EXPTIME", 1.0))
    saturation = float(header.get("DATASAT", 0.0))
    # Historical SECCHI Level-0 headers often encode DATASAT=0 as "not set".
    # Treating that sentinel as a threshold would mask the complete detector.
    bad = (data >= saturation) if saturation > 0 else np.zeros(data.shape, dtype=bool)
    data = (data - bias) / exptime
    data[bad] = np.nan
    return data, header


def world_to_pixel(header, radii, pa_deg):
    pa = np.deg2rad(pa_deg)[None, :]
    rsun = float(header["RSUN"])
    xw = radii[:, None] * rsun * np.sin(pa)
    yw = radii[:, None] * rsun * np.cos(pa)
    pc = np.array([[float(header.get("PC1_1", 1.0)), float(header.get("PC1_2", 0.0))],
                   [float(header.get("PC2_1", 0.0)), float(header.get("PC2_2", 1.0))]])
    crval = np.array([float(header.get("CRVAL1", 0.0)), float(header.get("CRVAL2", 0.0))])
    cdelt = np.array([float(header["CDELT1"]), float(header["CDELT2"])])
    delta = np.stack([xw - crval[0], yw - crval[1]], axis=0).reshape(2, -1)
    pix_delta = np.linalg.inv(pc) @ delta
    pix_delta /= cdelt[:, None]
    xx = float(header["CRPIX1"]) - 1.0 + pix_delta[0].reshape(xw.shape)
    yy = float(header["CRPIX2"]) - 1.0 + pix_delta[1].reshape(yw.shape)

    # Tangent vector d(x_world,y_world)/d(PA), transformed into detector pixels.
    tx = np.cos(pa)
    ty = -np.sin(pa)
    tangent = np.stack([np.broadcast_to(tx, xw.shape), np.broadcast_to(ty, yw.shape)], axis=0).reshape(2, -1)
    tangent_pix = np.linalg.inv(pc) @ tangent
    tangent_pix /= cdelt[:, None]
    chi = np.arctan2(tangent_pix[1], tangent_pix[0]).reshape(xw.shape)
    return xx, yy, chi


def process_triplet(commanded, radii, pa):
    date = commanded.strftime("%Y%m%d")
    hm = commanded.strftime("%H%M")
    cache_path = TRIPLET_CACHE / f"{date}_{hm}_fitpol_sector.npz"
    if cache_path.exists():
        try:
            with np.load(cache_path) as saved:
                observed = core.parse_utc(str(saved["observed_iso"].item()))
                return commanded, observed, saved["total_b"], saved["pb"], \
                    float(saved["fit_quality"].item()), str(saved["filenames"].item())
        except Exception:
            # A process interruption can leave only a partial ZIP container.
            cache_path.unlink(missing_ok=True)
    base = f"https://secchi.nrl.navy.mil/postflight/cor1/L0/a/seq/{date}/"
    names = [f"{date}_{hm}{sec}_s4c1A.fts" for sec in ("00", "09", "18")]
    blobs = [fetch(base + name) for name in names]
    images, headers = zip(*(read_image(blob) for blob in blobs))
    xx, yy, chi = world_to_pixel(headers[0], radii, pa)
    samples = np.stack([
        map_coordinates(image, [yy, xx], order=1, mode="constant", cval=np.nan)
        for image in images
    ], axis=0)
    polar_angles = np.deg2rad([float(h["POLAR"]) for h in headers])[:, None, None]
    c = np.cos(2.0 * (polar_angles - chi[None]))
    total_b = (2.0 / 3.0) * np.nansum(samples, axis=0)
    signed_pb = (4.0 / 3.0) * np.nansum(samples * c, axis=0)
    predicted = 0.5 * (total_b[None] + signed_pb[None] * c)
    residual = np.sqrt(np.nanmean((samples - predicted) ** 2, axis=0))
    scale = np.nanmedian(np.abs(samples), axis=0)
    fractional_residual = residual / np.maximum(scale, 1e-6)
    observed = core.parse_utc(str(headers[0]["DATE-OBS"]) + "Z")
    total_b = total_b.astype(np.float32)
    pb = np.abs(signed_pb).astype(np.float32)
    fit_quality = float(np.nanmedian(fractional_residual))
    filenames = "|".join(names)
    temporary = cache_path.with_name(cache_path.stem + ".tmp.npz")
    np.savez_compressed(
        temporary,
        observed_iso=np.array(observed.isoformat()),
        total_b=total_b,
        pb=pb,
        fit_quality=np.array(fit_quality),
        filenames=np.array(filenames),
    )
    temporary.replace(cache_path)
    return commanded, observed, total_b, pb, fit_quality, filenames


def available_sequence_times():
    """Resolve real complete triplets, then sample them near the frozen 15-min grid."""
    dates = []
    day = core.START.date()
    while day <= core.STOP.date():
        dates.append(day)
        day += timedelta(days=1)
    complete = []
    for day in dates:
        date = day.strftime("%Y%m%d")
        url = f"https://secchi.nrl.navy.mil/postflight/cor1/L0/a/seq/{date}/"
        html = fetch(url).decode("utf-8", "replace")
        names = set(re.findall(r'href="([^"]+_s4c1A\.fts)"', html))
        for name in sorted(names):
            match = re.match(r"(\d{8})_(\d{4})00_s4c1A\.fts$", name)
            if not match:
                continue
            prefix = f"{match.group(1)}_{match.group(2)}"
            if f"{prefix}09_s4c1A.fts" not in names or f"{prefix}18_s4c1A.fts" not in names:
                continue
            dt = core.parse_utc(
                f"{match.group(1)[:4]}-{match.group(1)[4:6]}-{match.group(1)[6:]}T"
                f"{match.group(2)[:2]}:{match.group(2)[2:]}:00Z"
            )
            if core.START <= dt <= core.STOP:
                complete.append(dt)
    complete.sort()
    targets = core.requested_times()
    chosen = []
    for target in targets:
        if not complete:
            break
        nearest = min(complete, key=lambda dt: abs((dt - target).total_seconds()))
        if abs((nearest - target).total_seconds()) <= 7.6 * 60 and nearest not in chosen:
            chosen.append(nearest)
    return chosen


def build_or_load_level0():
    if CACHE.exists():
        d = np.load(CACHE)
        return {k: d[k] for k in d.files}
    radii = np.arange(1.45, 3.351, 0.025, dtype=float)
    pa = np.arange(140.0, 210.001, 1.0, dtype=float)
    requests = available_sequence_times()
    print(f"Downloading and fitting {len(requests)} COR1 Level-0 triplets...", flush=True)
    results, failed = [], []
    with ThreadPoolExecutor(max_workers=16) as pool:
        futures = {pool.submit(process_triplet, dt, radii, pa): dt for dt in requests}
        for number, future in enumerate(as_completed(futures), 1):
            try:
                results.append(future.result())
            except Exception as exc:
                failed.append((futures[future].isoformat(), str(exc)))
                if len(failed) <= 5:
                    print(f"    failed {failed[-1][0]}: {failed[-1][1]}", flush=True)
            if number % 18 == 0 or number == len(requests):
                print(f"  Level-0 triplets: {number}/{len(requests)}; failed={len(failed)}", flush=True)
    results.sort(key=lambda item: item[1])
    if len(results) < 0.90 * len(requests):
        raise RuntimeError(f"Too many failed triplets: {len(failed)}")
    observed = [item[1] for item in results]
    minutes = np.array([(dt - core.T0).total_seconds() / 60.0 for dt in observed])
    total_b = np.stack([item[2] for item in results])
    pb = np.stack([item[3] for item in results])
    fit_quality = np.array([item[4] for item in results])
    filenames = np.array([item[5] for item in results])
    np.savez_compressed(
        CACHE, total_b=total_b, pb=pb, minutes=minutes, radii=radii, pa=pa,
        fit_quality=fit_quality, filenames=filenames,
        failed=np.array([json.dumps(x) for x in failed]),
        level=np.array("Level-0 bias/exposure normalized tangential two-parameter fitpol diagnostic"),
    )
    return {"total_b": total_b, "pb": pb, "minutes": minutes, "radii": radii,
            "pa": pa, "fit_quality": fit_quality, "filenames": filenames}


def regularize_good(data):
    quality = np.asarray(data["fit_quality"], float)
    center = np.nanmedian(quality)
    mad = 1.4826 * np.nanmedian(np.abs(quality - center))
    good = np.isfinite(quality) & (quality <= center + 6 * max(mad, 1e-8))
    grid = np.arange(math.ceil(np.min(data["minutes"][good]) / core.CADENCE_MIN) * core.CADENCE_MIN,
                     math.floor(np.max(data["minutes"][good]) / core.CADENCE_MIN) * core.CADENCE_MIN + 0.1,
                     core.CADENCE_MIN)
    variants = {}
    for name in ("total_b", "pb"):
        cube = np.asarray(data[name], float)
        flat = cube.reshape(len(cube), -1)
        out = np.empty((len(grid), flat.shape[1]), float)
        for j in range(flat.shape[1]):
            valid = good & np.isfinite(flat[:, j])
            out[:, j] = np.interp(grid, data["minutes"][valid], flat[valid, j]) if np.sum(valid) >= 2 else np.nan
        variants[name] = out.reshape((len(grid),) + cube.shape[1:])
    return grid, variants, good


def cross_validated_variants_at_shift(grid, zmaps, shifted_events, lags, families):
    even = np.arange(len(shifted_events)) % 2 == 0
    folds = ((even, ~even), (~even, even))
    output = {}
    for family_name, models in families.items():
        held_out, chosen = [], []
        for train, valid in folds:
            train_scores = []
            for product, zmap in zmaps.items():
                stack = core.event_stack(grid, zmap, shifted_events[train], lags)
                scores = core.sample_models(stack, lags, models)
                train_scores.append(scores)
            train_scores = np.stack(train_scores)
            product_names = list(zmaps)
            product_index, model_index = np.unravel_index(np.nanargmax(train_scores), train_scores.shape)
            product = product_names[int(product_index)]
            valid_stack = core.event_stack(grid, zmaps[product], shifted_events[valid], lags)
            held_out.append(core.score_fixed_model(valid_stack, lags, models[int(model_index)]))
            chosen.append((product, int(model_index)))
        output[family_name] = {"cv_score": float(np.mean(held_out)), "fold_scores": held_out,
                               "chosen": chosen}
    return output


def empirical_variant_tests(grid, zmaps, lags, families):
    observed = cross_validated_variants_at_shift(grid, zmaps, core.EVENT_TIMES, lags, families)
    min_lag = min(float(np.min(m.tau)) for group in families.values() for m in group)
    shift_lo = grid[0] - np.min(core.EVENT_TIMES) - min_lag
    shift_hi = grid[-1] - np.max(core.EVENT_TIMES)
    shifts = np.arange(math.ceil(shift_lo / core.CADENCE_MIN) * core.CADENCE_MIN,
                       math.floor(shift_hi / core.CADENCE_MIN) * core.CADENCE_MIN + 0.1,
                       core.CADENCE_MIN)
    null_shifts = shifts[np.abs(shifts) > 1e-6]
    null = {name: [] for name in families}
    gains = {"acceleration_minus_constant": [], "deceleration_minus_constant": []}
    for number, shift in enumerate(null_shifts, 1):
        result = cross_validated_variants_at_shift(grid, zmaps, core.EVENT_TIMES + shift, lags, families)
        for name in families:
            null[name].append(result[name]["cv_score"])
        gains["acceleration_minus_constant"].append(result["acceleration"]["cv_score"] - result["constant"]["cv_score"])
        gains["deceleration_minus_constant"].append(result["deceleration"]["cv_score"] - result["constant"]["cv_score"])
        if number % 30 == 0:
            print(f"  Level-0 null shifts: {number}/{len(null_shifts)}", flush=True)
    for name, models in families.items():
        values = np.asarray(null[name])
        obs = observed[name]["cv_score"]
        observed[name]["raw_p_common_shift"] = float((np.sum(values >= obs) + 1) / (len(values) + 1))
        observed[name]["null_mean"] = float(np.mean(values))
        observed[name]["null_std"] = float(np.std(values))
        observed[name]["null_scores"] = values.tolist()
        observed[name]["chosen_models"] = [
            {"product": product, **core.model_to_dict(models[index])}
            for product, index in observed[name]["chosen"]
        ]
    for key, values in gains.items():
        family = key.split("_minus_")[0]
        obs = observed[family]["cv_score"] - observed["constant"]["cv_score"]
        values = np.asarray(values)
        observed[key] = {"observed_score_gain": float(obs),
                         "raw_p_common_shift": float((np.sum(values >= obs) + 1) / (len(values) + 1)),
                         "null_mean": float(np.mean(values))}
    observed["null_shift_count"] = int(len(null_shifts))
    observed["null_shift_range_min"] = [float(np.min(null_shifts)), float(np.max(null_shifts))]
    return observed


def best_all_data(grid, zmaps, lags, families):
    output = {}
    for family_name, models in families.items():
        candidates = []
        for product, zmap in zmaps.items():
            stack = core.event_stack(grid, zmap, core.EVENT_TIMES, lags)
            scores = core.sample_models(stack, lags, models)
            i = int(np.nanargmax(scores))
            candidates.append((float(scores[i]), product, models[i], stack))
        score, product, model, stack = max(candidates, key=lambda x: x[0])
        output[family_name] = {"score": score, "product": product, "model": model, "stack": stack}
    return output


def make_figure(grid, cubes, radii_all, pa, path, zmaps, lags, best, tests):
    fig, axes = plt.subplots(2, 3, figsize=(16, 10), constrained_layout=True)
    static = np.nanmedian(cubes["pb"], axis=0)
    vmax = np.nanpercentile(static, 99)
    axes[0, 0].pcolormesh(pa, radii_all, static, shading="auto", cmap="gray", vmin=0, vmax=vmax)
    axes[0, 0].plot(path, radii_all, color="#e45756", lw=2)
    axes[0, 0].set(title="Level-0 tangential-fit pB and traced ray", xlabel="PA (deg)",
                   ylabel=r"Radius ($R_\odot$)")
    colors = {"constant": "white", "acceleration": "lime", "deceleration": "gold"}
    display = best["acceleration"] if tests["acceleration"]["raw_p_common_shift"] <= tests["deceleration"]["raw_p_common_shift"] else best["deceleration"]
    stack = display["stack"]
    lim = np.nanpercentile(np.abs(stack), 98)
    axes[0, 1].pcolormesh(lags, core.RADII_WORK, stack.T, shading="auto", cmap="RdBu_r", vmin=-lim, vmax=lim)
    for name, row in best.items():
        axes[0, 1].plot(row["model"].tau, core.RADII_WORK, color=colors[name], lw=2,
                        label=f"{name} ({row['product']})")
    axes[0, 1].set(title=f"12-event stack shown in {display['product']}",
                   xlabel=r"Minutes relative to 3 R$_\odot$ anchor", ylabel=r"Radius ($R_\odot$)")
    axes[0, 1].legend(fontsize=8, loc="lower right")
    names = ["constant", "acceleration", "deceleration"]
    pvals = [tests[n]["raw_p_common_shift"] for n in names]
    bars = axes[0, 2].bar(names, -np.log10(pvals), color=[colors[n] for n in names], edgecolor="k")
    axes[0, 2].axhline(-np.log10(0.05), color="k", ls="--")
    axes[0, 2].set(title="Held-out Level-0 test (no BH)", ylabel=r"$-\log_{10}p_{raw}$")
    for bar, p in zip(bars, pvals):
        axes[0, 2].text(bar.get_x()+bar.get_width()/2, bar.get_height()+.03, f"p={p:.3f}", ha="center")
    for name, row in best.items():
        model = row["model"]
        speed = np.where(core.RADII_WORK < (model.change_radius or -99), model.v_inner, model.v_outer)
        axes[1, 0].step(core.RADII_WORK, speed, where="mid", color=colors[name], lw=2,
                        label=f"{name}: {row['product']}")
    axes[1, 0].set(title="Best all-event profiles (visualization only)", xlabel=r"Radius ($R_\odot$)",
                   ylabel="Pattern speed (km s$^{-1}$)")
    axes[1, 0].legend(fontsize=8)
    chosen_family = min(names, key=lambda n: tests[n]["raw_p_common_shift"])
    chosen = best[chosen_family]
    amplitude, jitter, offsets = core.profile_diagnostics(grid, zmaps[chosen["product"]], core.RADII_WORK, chosen["model"])
    axes[1, 1].plot(core.RADII_WORK, amplitude, color="#2a9d8f", lw=2)
    if chosen["model"].change_radius:
        axes[1, 1].axvline(chosen["model"].change_radius, color=colors[chosen_family], ls="--")
    axes[1, 1].axhline(0, color="0.5", lw=1)
    axes[1, 1].set(title=f"Ridge amplitude: {chosen_family}, {chosen['product']}",
                   xlabel=r"Radius ($R_\odot$)", ylabel="Robust BFF z")
    axes[1, 2].plot(core.RADII_WORK, jitter, color="#6a4c93", lw=2)
    if chosen["model"].change_radius:
        axes[1, 2].axvline(chosen["model"].change_radius, color=colors[chosen_family], ls="--")
    axes[1, 2].set(title="Event-to-event timing jitter", xlabel=r"Radius ($R_\odot$)", ylabel="Robust scatter (min)")
    fig.suptitle("11–14 January 2008: Level-0 nonlinear COR1 transport test", fontsize=16)
    fig.savefig(OUT / "pds_20080111_14_level0_nonlinear_transport_diagnostic.png", dpi=180)
    plt.close(fig)
    return chosen_family, chosen, amplitude, jitter, offsets


def write_report(data, good, tests, best, chosen_family, chosen, amplitude, jitter):
    accel = tests["acceleration"]["chosen_models"]
    decel = tests["deceleration"]["chosen_models"]
    accel_stable = abs(accel[0]["change_radius_rsun"] - accel[1]["change_radius_rsun"]) <= .2
    decel_stable = abs(decel[0]["change_radius_rsun"] - decel[1]["change_radius_rsun"]) <= .2
    selected_stable = accel_stable if chosen_family == "acceleration" else decel_stable if chosen_family == "deceleration" else True
    gain_key = f"{chosen_family}_minus_constant" if chosen_family != "constant" else None
    candidate = (chosen_family != "constant" and tests[chosen_family]["raw_p_common_shift"] < .05
                 and tests[gain_key]["raw_p_common_shift"] < .05 and selected_stable)
    lines = [
        "# Level-0 nonlinear PDS transport test: 11–14 January 2008",
        "",
        "## Main result",
        "",
    ]
    if candidate:
        lines += [
            f"A cross-validated {chosen_family} change is detected in COR1 timing, but it is",
            "only a kinematic transition candidate.  A shock interpretation still requires",
            "absolute Level-1 pB compression and an independent M_f > 1 test.",
        ]
    else:
        lines += [
            "No stationary, common acceleration or deceleration break is established by the",
            "held-out COR1 ridges.  Allowing a nonlinear travel-time law does not rescue a",
            "single deterministic EUVI/COR1-to-COR2 propagation chain.",
        ]
    lines += [
        "",
        "## Frozen design",
        "",
        "- 12 z >= 1.25 COR2 anchors fixed independently at 3 R_sun.",
        f"- Fixed seed PA = {core.SEED_PA_DEG:.1f} deg.",
        "- COR1 model families fixed in advance: constant speed, one acceleration break,",
        "  or one deceleration break.",
        "- Each fold chooses tB or pB, speed(s), and change radius on six training ridges;",
        "  the choice is then frozen and scored on the other six ridges.",
        "- The entire selection is repeated in every common-shift null realization.",
        "- Raw empirical p values are reported; no BH adjustment.",
        "",
        "## Data",
        "",
        f"- Successfully fitted Level-0 triplets: {len(data['minutes'])}; quality-retained: {int(np.sum(good))}.",
        "- Each triplet is bias-subtracted and exposure-normalized.",
        "- tB and the tangential polarized component are fitted from the 0/120/240 degree exposures.",
        "- This diagnostic does not apply the complete SECCHI_PREP calibration/background chain;",
        "  absolute compression ratios must not be inferred from it.",
        "",
        "## Held-out results",
        "",
        "| Family | CV score | Raw common-shift p | Fold 1 choice | Fold 2 choice |",
        "| --- | ---: | ---: | --- | --- |",
    ]
    for name in ("constant", "acceleration", "deceleration"):
        choices = []
        for m in tests[name]["chosen_models"]:
            if m["change_radius_rsun"] is None:
                choices.append(f"{m['product']}, {m['v_inner_km_s']:.0f} km/s")
            else:
                choices.append(f"{m['product']}, {m['v_inner_km_s']:.0f}->{m['v_outer_km_s']:.0f} km/s at {m['change_radius_rsun']:.2f} R_sun")
        lines.append(f"| {name} | {tests[name]['cv_score']:.4f} | {tests[name]['raw_p_common_shift']:.4f} | {choices[0]} | {choices[1]} |")
    lines += [
        "",
        f"- Acceleration gain over constant: {tests['acceleration_minus_constant']['observed_score_gain']:.4f}; "
        f"raw p={tests['acceleration_minus_constant']['raw_p_common_shift']:.4f}.",
        f"- Deceleration gain over constant: {tests['deceleration_minus_constant']['observed_score_gain']:.4f}; "
        f"raw p={tests['deceleration_minus_constant']['raw_p_common_shift']:.4f}.",
        f"- Acceleration change radius stable within 0.2 R_sun across folds: {accel_stable}.",
        f"- Deceleration change radius stable within 0.2 R_sun across folds: {decel_stable}.",
        "",
        "## All-event fits (visualization, not significance)",
        "",
    ]
    for name, row in best.items():
        m = row["model"]
        if m.change_radius is None:
            detail = f"v={m.v_inner:.0f} km/s"
        else:
            detail = f"{m.v_inner:.0f}->{m.v_outer:.0f} km/s at {m.change_radius:.2f} R_sun"
        lines.append(f"- {name}: {row['product']}, {detail}, stacked score={row['score']:.4f}.")
    lines += [
        "",
        "## Physical interpretation",
        "",
        "The fixed 200 km/s test was too restrictive, but relaxing it is not sufficient by",
        "itself.  A stationary shock/nozzle transition must predict independent ridges with",
        "one reproducible travel-time law.  If that condition fails, the remaining physical",
        "possibilities are local release near the cusp, a moving/time-dependent throat or",
        "shock that produces event-dependent delays, and a height-dependent visibility effect.",
        "",
        "A true shock claim additionally requires a repeatable density jump at the same height,",
        "M_f=(u_n-V_sh,n)/c_f > 1 upstream, and Rankine–Hugoniot consistency.  None of those",
        "conditions is inferred from the timing p value alone.",
    ]
    (OUT / "PDS_20080111_14_Level0_nonlinear_transport_no_BH.md").write_text("\n".join(lines)+"\n", encoding="utf-8")
    return candidate, accel_stable, decel_stable


def main():
    data = build_or_load_level0()
    grid, cubes, good = regularize_good(data)
    path, static = core.trace_path(cubes["pb"], data["radii"], data["pa"], core.SEED_PA_DEG)
    maps = {}
    raw_maps = {}
    use = (data["radii"] >= 1.55) & (data["radii"] <= 3.001)
    core.RADII_WORK = data["radii"][use]
    for name, cube in cubes.items():
        values = core.sample_path(cube, data["pa"], path, half_width=2.5)
        raw, z = core.bff_maps(values)
        maps[name] = z[:, use]
        raw_maps[name] = raw[:, use]
    lags = np.arange(-920.0, 0.1, core.CADENCE_MIN)
    families = core.model_families(core.RADII_WORK)
    print("Level-0 model counts:", {k: len(v) for k, v in families.items()}, flush=True)
    tests = empirical_variant_tests(grid, maps, lags, families)
    best = best_all_data(grid, maps, lags, families)
    chosen_family, chosen, amplitude, jitter, offsets = make_figure(
        grid, cubes, data["radii"], data["pa"], path, maps, lags, best, tests
    )
    candidate, accel_stable, decel_stable = write_report(
        data, good, tests, best, chosen_family, chosen, amplitude, jitter
    )
    output = {
        "event": "2008-01-11--14",
        "data_level": "Level-0 bias/exposure normalized tangential fitpol diagnostic",
        "frames": int(len(data["minutes"])), "quality_retained": int(np.sum(good)),
        "cor2_fixed_pa_deg": core.SEED_PA_DEG,
        "strong_event_count": int(len(core.EVENT_TIMES)),
        "tests": tests,
        "best_all_event_models": {
            name: {"product": row["product"], **core.model_to_dict(row["model"]),
                   "stacked_score": row["score"]}
            for name, row in best.items()
        },
        "diagnostic_selected_family": chosen_family,
        "diagnostic_selected_product": chosen["product"],
        "diagnostic_selected_model": core.model_to_dict(chosen["model"]),
        "acceleration_change_radius_stable": accel_stable,
        "deceleration_change_radius_stable": decel_stable,
        "strict_stationary_transition_candidate": candidate,
        "radii_rsun": core.RADII_WORK.tolist(),
        "ridge_amplitude_profile": amplitude.tolist(),
        "timing_jitter_profile_min": jitter.tolist(),
        "event_timing_offsets_min": offsets.tolist(),
        "traced_path_pa_deg": path.tolist(),
        "no_BH": True,
        "calibration_caveat": "Not full SECCHI_PREP; no absolute pB compression ratio",
    }
    (OUT / "pds_20080111_14_level0_nonlinear_transport_results.json").write_text(
        json.dumps(output, indent=2), encoding="utf-8"
    )
    np.savez_compressed(
        OUT / "pds_20080111_14_level0_nonlinear_transport_maps.npz",
        minutes=grid, radii=core.RADII_WORK, pa=data["pa"], traced_path_pa=path,
        z_total_b=maps["total_b"].astype(np.float32), z_pb=maps["pb"].astype(np.float32),
        raw_total_b=raw_maps["total_b"].astype(np.float32), raw_pb=raw_maps["pb"].astype(np.float32),
        event_times=core.EVENT_TIMES, lags=lags,
    )
    print(json.dumps({
        "p": {n: tests[n]["raw_p_common_shift"] for n in ("constant", "acceleration", "deceleration")},
        "gain_p": {"acceleration": tests["acceleration_minus_constant"]["raw_p_common_shift"],
                   "deceleration": tests["deceleration_minus_constant"]["raw_p_common_shift"]},
        "fold_models": {n: tests[n]["chosen_models"] for n in ("constant", "acceleration", "deceleration")},
        "strict_stationary_transition_candidate": candidate,
    }, indent=2))


if __name__ == "__main__":
    main()
