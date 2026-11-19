#!/usr/bin/env python3
"""Build vector figures and compact LaTeX tables for the ApJ manuscript.

The script reads only frozen, tabular or compact NPZ analysis products.  It
does not recalculate event statistics or inspect the raw images.  All figure
panels are written as vector PDF files.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, Ellipse, FancyArrowPatch, Polygon, Rectangle


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
FIG = ROOT / "figures"
TABLES = ROOT / "tables"

COLORS = {
    "navy": "#17365D",
    "blue": "#2F6DA1",
    "cyan": "#5AA6A6",
    "orange": "#D17A22",
    "red": "#A43A3A",
    "green": "#497A55",
    "gray": "#6A7075",
    "light": "#E8EDF2",
    "ink": "#20252A",
}


def style() -> None:
    mpl.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 8.2,
        "axes.labelsize": 8.2,
        "axes.titlesize": 9.0,
        "xtick.labelsize": 7.3,
        "ytick.labelsize": 7.3,
        "legend.fontsize": 7.2,
        "axes.linewidth": 0.7,
        "xtick.major.width": 0.7,
        "ytick.major.width": 0.7,
        "savefig.bbox": "tight",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })


def read_csv(name: str) -> list[dict[str, str]]:
    with (DATA / name).open(newline="") as handle:
        return list(csv.DictReader(handle))


def panel_label(ax, label: str) -> None:
    ax.text(0.01, 0.98, label, transform=ax.transAxes, va="top", ha="left",
            fontsize=9.5, fontweight="bold")


def save(fig, name: str) -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / name, format="pdf", dpi=300)
    plt.close(fig)


def figure1_concept() -> None:
    fig, axs = plt.subplots(2, 2, figsize=(7.15, 5.35), constrained_layout=True)
    snapshot_titles = [r"$t=t_0$", r"$t=t_0+P$", r"$t=t_0+2P$", r"$t=t_0+3P$"]
    for n, (ax, title) in enumerate(zip(axs.flat, snapshot_titles)):
        ax.set_aspect("equal")
        ax.set_xlim(-0.25, 7.2)
        ax.set_ylim(-1.55, 1.55)
        ax.axis("off")
        ax.set_title(title + r",  $P=120$ min (illustrative)", loc="left",
                     fontsize=8.2, weight="bold", color=COLORS["navy"])

        # Solar limb and compact closed arcade.
        ax.add_patch(Circle((0.20, 0), 0.72, facecolor="#F6C85F",
                            edgecolor=COLORS["orange"], lw=0.9))
        xloop = np.linspace(0.70, 2.42, 120)
        for height in (0.42, 0.68, 0.94):
            yloop = height * np.sin(np.pi * (xloop - 0.70) / 1.72)
            ax.plot(xloop, yloop, color=COLORS["gray"], lw=0.65)
            ax.plot(xloop, -yloop, color=COLORS["gray"], lw=0.65)
        ax.add_patch(Ellipse((1.37, 0), 1.22, 0.55, fill=False,
                             edgecolor=COLORS["blue"], lw=1.25, ls="--"))
        ax.text(1.37, 1.20, "compact MHD\nresonator candidate", ha="center", va="top",
                fontsize=6.7, color=COLORS["blue"], weight="bold")

        # Cusp/current-sheet gate and expanding guide.
        cusp_x = 2.42
        ax.plot([cusp_x, 3.08], [0, 0], color=COLORS["red"], lw=1.0, ls="--")
        ax.scatter([cusp_x], [0], s=18, color=COLORS["red"], zorder=5)
        ax.text(2.65, -1.05, "reforming\ngate", ha="center", fontsize=6.7,
                color=COLORS["red"], weight="bold")
        xs = np.linspace(cusp_x, 7.0, 160)
        width = 0.15 + 0.115 * (xs - cusp_x) ** 0.82
        ax.plot(xs, width, color=COLORS["navy"], lw=1.05)
        ax.plot(xs, -width, color=COLORS["navy"], lw=1.05)
        ax.fill_between(xs, -width, width, color=COLORS["blue"], alpha=0.06)
        ax.add_patch(FancyArrowPatch((2.75, 0), (6.85, 0), arrowstyle="-|>",
                                     mutation_scale=9, color=COLORS["ink"], lw=0.65))

        # A newly released packet appears every cycle; older packets move outward
        # and broaden. This illustrates cadence-to-spacing conversion.
        for release in range(n + 1):
            age = n - release
            xp = 2.90 + 0.98 * age + 0.18 * age ** 2
            wp = 0.10 + 0.045 * age
            hp = 0.32 + 0.07 * age
            poly = Polygon([[xp - wp, -0.18], [xp, -hp], [xp + wp, -0.18],
                            [xp + wp, 0.18], [xp, hp], [xp - wp, 0.18]],
                           closed=True, facecolor=COLORS["cyan"],
                           edgecolor=COLORS["blue"], lw=0.75, alpha=0.85)
            ax.add_patch(poly)
        ax.annotate("new release", xy=(2.90, 0.38), xytext=(3.45, 0.96),
                    fontsize=6.5, ha="center", color=COLORS["green"],
                    arrowprops=dict(arrowstyle="->", lw=0.6, color=COLORS["green"]))
        if n == 3:
            ax.text(5.35, -1.22, "advected PDS train", ha="center",
                    fontsize=6.8, color=COLORS["navy"], weight="bold")

        panel_label(ax, f"({chr(97+n)})")
    fig.suptitle("Hypothesis: a 120-minute source modulation can release a spatially expanding PDS train",
                 fontsize=9.2, weight="bold", color=COLORS["ink"])
    save(fig, "fig1_observational_logic.pdf")


def figure6_synthesis() -> None:
    """Three transfer cartoons and the inference supported by the 12-event sample."""
    fig, axs = plt.subplots(1, 3, figsize=(7.15, 3.25), constrained_layout=True)
    t = np.linspace(0.8, 9.2, 120)
    colors = [COLORS["green"], COLORS["red"], COLORS["orange"]]
    titles = ["Smooth stationary nozzle", "Stationary shock", "Reforming MHD gate"]
    for k, (ax, color, title) in enumerate(zip(axs, colors, titles)):
        ax.set_xlim(0, 10)
        ax.set_ylim(1.45, 3.05)
        ax.set_xlabel("time")
        if k == 0:
            ax.set_ylabel(r"height ($R_\odot$)")
        else:
            ax.set_yticklabels([])
        ax.set_title(title, color=color, weight="bold", fontsize=8.3)
        ax.spines[["top", "right"]].set_visible(False)
        for i in range(4):
            phase = i * 1.75
            if k < 2:
                r = 1.55 + 0.145 * (t - phase) + 0.010 * (t - phase) ** 2
            else:
                jitter = [0.00, 0.18, -0.12, 0.26][i]
                r = 1.55 + 0.145 * (t - phase - jitter) + 0.010 * (t - phase - jitter) ** 2
                r += 0.035 * np.sin(1.7 * t + i)
            keep = (r >= 1.5) & (r <= 3.0)
            ax.plot(t[keep], r[keep], lw=1.15, color=COLORS["navy"], alpha=0.92)
        if k == 0:
            ax.text(5.1, 2.88, "parallel tracks\ngradual broadening", ha="center", va="top", fontsize=7.0)
        elif k == 1:
            ax.axhline(2.35, color=COLORS["red"], lw=1.2, ls="--")
            ax.text(5.0, 2.92, "fixed-height jump;\ncadence may survive", ha="center", va="top", fontsize=7.0)
        else:
            gate = 2.35 + 0.12 * np.sin(0.75 * t)
            ax.plot(t, gate, color=COLORS["orange"], lw=1.2, ls="--")
            ax.text(5.0, 2.92, "order survives;\ndelays differ", ha="center", va="top", fontsize=7.0)
        ax.tick_params(labelsize=6.8)
    fig.text(0.5, -0.01,
             "Observed combination: high-COR1 localization + expansion-stable order + no closed common-delay or jump test",
             ha="center", va="top", fontsize=7.7, color=COLORS["ink"])
    save(fig, "fig6_physical_synthesis.pdf")


def figure2_survey() -> None:
    rows = read_csv("all12_events.csv")
    ev = np.arange(1, 13)
    peaks = np.array([int(r["strict_matched_peak_count"]) for r in rows])
    pcell = np.array([float(r["strict_cell_stat_p_raw"]) for r in rows])
    pbest = np.array([float(r["best_ridge_p_raw"]) for r in rows])
    pcommon = np.array([float(r["common_ridge_p_raw"]) for r in rows])
    phase = np.array([float(r["tb_pb_global_phase_difference_min"]) for r in rows])
    cat = [r["category"] for r in rows]
    cmap = {"A": COLORS["green"], "B": COLORS["blue"], "C": COLORS["orange"], "D": COLORS["gray"]}
    colors = [cmap[c] for c in cat]

    fig, axs = plt.subplots(2, 2, figsize=(7.15, 5.45), constrained_layout=True)
    ax = axs[0, 0]
    ax.bar(ev, peaks, color=colors, edgecolor="white", linewidth=0.5)
    ax.axhline(3, color=COLORS["red"], ls="--", lw=0.8)
    ax.set(ylabel="strict matched peaks", xlabel="event number", xticks=ev, ylim=(0, 4.7))
    panel_label(ax, "(a)")

    ax = axs[0, 1]
    floor = 1e-3
    ax.plot(ev, -np.log10(np.maximum(pcell, floor)), "o-", label="cell statistic", color=COLORS["navy"], ms=3.5)
    ax.plot(ev, -np.log10(np.maximum(pbest, floor)), "s--", label="best ridge", color=COLORS["orange"], ms=3.2)
    ax.plot(ev, -np.log10(np.maximum(pcommon, floor)), "^:", label="common ridge", color=COLORS["cyan"], ms=3.4)
    ax.axhline(-np.log10(0.05), color=COLORS["red"], ls="--", lw=0.8, label="raw p=0.05")
    ax.set(ylabel=r"$-\log_{10}(p_{\rm raw})$", xlabel="event number", xticks=ev, ylim=(0, 1.75))
    ax.legend(ncol=2, frameon=False, loc="upper left")
    panel_label(ax, "(b)")

    ax = axs[1, 0]
    ax.bar(ev, phase, color=colors, edgecolor="white", linewidth=0.5)
    ax.axhspan(0, 30, color=COLORS["green"], alpha=0.10)
    ax.axhline(30, color=COLORS["green"], lw=0.8, ls="--")
    ax.set(ylabel=r"global $tB$--$pB$ offset (min)", xlabel="event number", xticks=ev, ylim=(0, 195))
    panel_label(ax, "(c)")

    ax = axs[1, 1]
    counts = [cat.count(c) for c in "ABCD"]
    ax.bar(np.arange(4), counts, color=[cmap[c] for c in "ABCD"], width=0.72)
    for x, n in enumerate(counts):
        ax.text(x, n+0.10, f"{n}/12", ha="center", va="bottom", fontsize=8)
    labels = ["A\nordered", "B\ncompatible", "C\nweak", "D\nincomplete"]
    ax.set(xticks=np.arange(4), xticklabels=labels, ylabel="events", ylim=(0, 5.1))
    panel_label(ax, "(d)")
    save(fig, "fig2_event_survey.pdf")


def robust_limits(a: np.ndarray) -> tuple[float, float]:
    q1, q99 = np.nanpercentile(a, [2, 98])
    m = max(abs(q1), abs(q99), 1e-6)
    return -m, m


def figure3_maps() -> None:
    z = np.load(DATA / "two_event_maps.npz")
    radii = z["radii"]
    pa = z["pa_offsets"]
    geo = json.loads((DATA / "spatial_expanding_results.json").read_text())
    products = [("event9_base60_total_b", "Event 9, $tB$"),
                ("event9_base60_pb", "Event 9, $pB$"),
                ("event12_base60_total_b", "Event 12, $tB$"),
                ("event12_base60_pb", "Event 12, $pB$")]
    fig, axs = plt.subplots(2, 2, figsize=(7.15, 5.45), sharex=True, sharey=True, constrained_layout=True)
    cmap = mpl.colormaps["RdBu_r"]
    for idx, (ax, (key, title)) in enumerate(zip(axs.flat, products)):
        arr = np.asarray(z[key], float)
        lo, hi = robust_limits(arr)
        ax.pcolormesh(pa, radii, arr, cmap=cmap, shading="auto", vmin=lo, vmax=hi,
                      rasterized=False)
        event_key = "event9" if "event9" in key else "event12"
        gs = geo["events"]["measured"][event_key]["base60"]["geometries"]
        q = np.linspace(pa.min(), pa.max(), 120)
        for j, g in enumerate(gs):
            rv = float(g["vertex_radius_rsun"])
            sl = float(g["absolute_slope_rsun_per_degree"])
            col = COLORS["ink"] if j % 2 == 0 else COLORS["green"]
            ax.plot(q, rv + sl*q, color=col, lw=0.75)
            ax.plot(q, rv - sl*q, color=col, lw=0.75)
            ax.plot(0, rv, marker="o", ms=2.5, color=col)
        ax.axvline(0, color="white", lw=0.55, alpha=0.8)
        ax.set_title(title)
        ax.set_xlim(pa.min(), pa.max())
        ax.set_ylim(radii.min(), radii.max())
        # A patch-based color key remains vector in the output PDF.
        cax = ax.inset_axes([1.02, 0.0, 0.035, 1.0])
        bins = 18
        for k in range(bins):
            cax.add_patch(Rectangle((0, k), 1, 1, facecolor=cmap((k + 0.5) / bins),
                                    edgecolor="none"))
        cax.set(xlim=(0, 1), ylim=(0, bins), xticks=[],
                yticks=[0, bins/2, bins], yticklabels=[f"{lo:.1f}", "0", f"{hi:.1f}"])
        cax.yaxis.tick_right()
        cax.yaxis.set_label_position("right")
        cax.set_ylabel("robust normalized contrast", labelpad=3)
        cax.tick_params(length=2, pad=1)
        panel_label(ax, f"({chr(97+idx)})")
    for ax in axs[:, 0]:
        ax.set_ylabel(r"radius ($R_\odot$)")
    for ax in axs[1, :]:
        ax.set_xlabel("position-angle offset (deg)")
    save(fig, "fig3_candidate_maps.pdf")


def figure4_expansion() -> None:
    d = json.loads((DATA / "spatial_expanding_results.json").read_text())
    reps = ["bff", "base60", "nrgf60"]
    schemes = ["fixed", "spherical", "measured"]
    markers = {"bff": "o", "base60": "s", "nrgf60": "^"}
    rep_labels = {"bff": "BFF", "base60": "base-60", "nrgf60": "NRGF-60"}
    fig, axs = plt.subplots(1, 2, figsize=(7.15, 3.45), constrained_layout=True)

    ax = axs[0]
    xpos = np.arange(len(schemes))
    for event, color, shift in [("event9", COLORS["blue"], -0.07), ("event12", COLORS["orange"], 0.07)]:
        for rep in reps:
            yorder = [d["events"][s][event][rep]["p_outward_order"] for s in schemes]
            ytrans = [d["events"][s][event][rep]["p_transport_25_300"] for s in schemes]
            ax.plot(xpos + shift, yorder, marker=markers[rep], ms=4, lw=0.8, color=color, alpha=0.88)
            ax.plot(xpos + shift, ytrans, marker=markers[rep], ms=4, lw=0.8, color=color,
                    alpha=0.88, ls="--", fillstyle="none")
    ax.axhline(0.5, color=COLORS["red"], ls="--", lw=0.8)
    ax.set(xticks=xpos, xticklabels=["fixed", r"$r^2$", "measured"], ylim=(0, 1.02),
           ylabel="probability", xlabel="ridge-width law")
    ax.text(0.03, 0.08, "solid: outward order\ndashed: transport in 25--300 km s$^{-1}$",
            transform=ax.transAxes, fontsize=7.1,
            bbox=dict(fc="white", ec=COLORS["light"], alpha=0.9, pad=2))
    ax.text(0.69, 0.94, "blue: event 9\norange: event 12", transform=ax.transAxes,
            ha="left", va="top", fontsize=7.1)
    panel_label(ax, "(a)")

    ax = axs[1]
    pair = d["pair_null"]
    for i, scheme in enumerate(schemes):
        vals = [pair[scheme]["representations"][rep]["p_raw"] for rep in reps]
        primary = pair[scheme]["primary_pair_p_raw"]
        for j, (rep, val) in enumerate(zip(reps, vals)):
            ax.scatter(i + (j-1)*0.10, val, marker=markers[rep], s=30,
                       color=[COLORS["navy"], COLORS["cyan"], COLORS["green"]][j],
                       label=rep_labels[rep] if i == 0 else None, zorder=3)
        ax.scatter(i, primary, marker="*", s=75, color=COLORS["red"],
        label=r"primary $\geq$2/3" if i == 0 else None, zorder=4)
    ax.axhline(0.05, color=COLORS["red"], ls="--", lw=0.8)
    ax.set(xticks=xpos, xticklabels=["fixed", r"$r^2$", "measured"],
           ylabel="shifted-pair raw $p$", xlabel="ridge-width law", ylim=(0, 0.20))
    ax.legend(frameon=False, ncol=2, loc="upper left")
    panel_label(ax, "(b)")
    save(fig, "fig4_expansion_sensitivity.pdf")


def figure5_height_controls() -> None:
    rows = read_csv("height_matrix_no_BH.csv")
    order = ["EUVI 1.02–1.25", "EUVI 1.25–1.45", "COR1 1.55–2.00",
             "COR1 2.00–2.50", "COR1 2.50–3.00"]
    cases = ["2008-01-11–14", "2008-01-18–20", "2008-05-29–31", "2008-07-26–28"]
    case_labels = ["Jan 11--14", "Jan 18--20", "May 29--31", "Jul 26--28"]
    colors = [COLORS["navy"], COLORS["orange"], COLORS["gray"], COLORS["green"]]
    marks = ["o", "s", "^", "D"]
    fig, ax = plt.subplots(figsize=(7.15, 3.55))
    for case, label, color, marker in zip(cases, case_labels, colors, marks):
        sub = {r["height_band"]: float(r["raw_p_common_shift"]) for r in rows if r["case"] == case}
        vals = [sub[o] for o in order]
        ax.semilogy(np.arange(5), vals, marker=marker, ms=4.4, lw=1.05, color=color, label=label)
    ax.axhline(0.05, color=COLORS["red"], ls="--", lw=0.8, label="raw p=0.05")
    ax.axvspan(1.5, 4.5, color=COLORS["blue"], alpha=0.045)
    ax.set(xticks=np.arange(5), xticklabels=["EUVI\n1.02--1.25", "EUVI\n1.25--1.45",
                                            "COR1\n1.55--2.00", "COR1\n2.00--2.50", "COR1\n2.50--3.00"],
           ylabel="common-shift raw $p$", xlabel=r"height band ($R_\odot$)", ylim=(0.01, 1.15))
    ax.legend(frameon=False, ncol=3, loc="upper center")
    ax.text(4, 0.02479, "  Jan 11--14: 0.0248", va="center", color=COLORS["navy"], fontsize=7.4)
    panel_label(ax, "(a)")
    save(fig, "fig5_height_controls.pdf")


def figure7_resonator_height_chain() -> None:
    """Place the local source candidate and failed/positive links on one page."""
    rows = read_csv("resonator_candidate_metrics.csv")
    fig, axs = plt.subplots(1, 2, figsize=(7.15, 4.25), constrained_layout=True,
                            gridspec_kw={"width_ratios": [0.92, 1.35]})

    # Left: one vertical physical scale for the three SECCHI domains.
    ax = axs[0]
    ax.set_xlim(0, 1)
    ax.set_ylim(1.0, 3.25)
    ax.set_xticks([])
    ax.set_ylabel(r"heliocentric height ($R_\odot$)")
    ax.spines[["top", "right", "bottom"]].set_visible(False)
    zones = [
        (1.02, 1.45, "EUVI", COLORS["green"]),
        (1.55, 3.00, "COR1", COLORS["orange"]),
        (3.00, 3.22, "COR2 anchors", COLORS["blue"]),
    ]
    for y0, y1, label, color in zones:
        ax.add_patch(Rectangle((0.06, y0), 0.24, y1-y0, facecolor=color,
                               edgecolor="none", alpha=0.18))
        if label.startswith("COR2"):
            ax.text(0.18, (y0+y1)/2, "COR2", ha="center", va="center",
                    color=color, weight="bold", fontsize=6.7)
        else:
            ax.text(0.18, (y0+y1)/2, label, rotation=90, ha="center", va="center",
                    color=color, weight="bold", fontsize=7.2)

    ax.add_patch(Rectangle((0.38, 1.10), 0.52, 0.10, facecolor=COLORS["green"],
                           edgecolor=COLORS["green"], alpha=0.78))
    ax.text(0.64, 1.25, "EUVI 171 Å\n80--130 min excess\n$p_{raw}=0.0279$",
            ha="center", va="bottom", fontsize=7.0, color=COLORS["ink"])

    ax.add_patch(Rectangle((0.38, 2.50), 0.52, 0.50, facecolor=COLORS["orange"],
                           edgecolor=COLORS["orange"], alpha=0.56))
    ax.text(0.64, 2.73, "upper-COR1\nensemble localization",
            ha="center", va="center", fontsize=7.0, color=COLORS["ink"])

    ax.scatter([0.64], [3.08], marker="v", s=38, color=COLORS["blue"], zorder=4)
    ax.text(0.64, 3.16, "12 COR2-selected events", ha="center", va="bottom",
            fontsize=7.0, color=COLORS["blue"], weight="bold")

    ax.add_patch(FancyArrowPatch((0.64, 1.31), (0.64, 2.42), arrowstyle="-|>",
                                 mutation_scale=9, lw=0.9, ls="--",
                                 color=COLORS["red"]))
    ax.text(0.69, 1.86, "phase bridge\nnot closed", ha="left", va="center",
            fontsize=6.7, color=COLORS["red"], weight="bold")
    panel_label(ax, "(a)")
    ax.set_title("Where the evidence appears", fontsize=8.6, weight="bold")

    # Right: distinct empirical questions shown on a common -log10(p) display.
    ax = axs[1]
    labels = [
        "EUVI 171\nlocal band",
        "EUVI→COR1\nfixed continuation",
        "COR1 strong-event\naggregate",
        "EUVI→COR2\nfixed continuation",
    ]
    pvals = np.array([float(r["p_raw"]) for r in rows])
    scores = -np.log10(np.clip(pvals, 1e-6, 1.0))
    colors = [COLORS["green"], COLORS["gray"], COLORS["orange"], COLORS["gray"]]
    y = np.arange(len(rows))
    ax.barh(y, scores, color=colors, alpha=0.86, height=0.62)
    ax.axvline(-np.log10(0.05), color=COLORS["red"], lw=0.85, ls="--")
    ax.text(-np.log10(0.05)+0.03, 3.48, "raw p=0.05", color=COLORS["red"],
            fontsize=6.7, va="top")
    for yi, (score, p) in enumerate(zip(scores, pvals)):
        ax.text(score + 0.04, yi, f"p={p:.4g}", va="center", fontsize=7.0)
    ax.set(yticks=y, yticklabels=labels, xlabel=r"$-\log_{10}(p_{\rm raw})$",
           xlim=(0, 2.55))
    ax.invert_yaxis()
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="x", color=COLORS["light"], lw=0.6)
    ax.set_axisbelow(True)
    ax.set_title("Which observational links close?", fontsize=8.6, weight="bold")
    ax.text(0.02, -0.22,
            "Different tests; no multiplication of p-values.\nThe positive source and COR1 results are separated by a failed phase link.",
            transform=ax.transAxes, ha="left", va="top", fontsize=6.8,
            color=COLORS["gray"])
    panel_label(ax, "(b)")
    save(fig, "fig7_resonator_height_chain.pdf")


def write_event_table() -> None:
    rows = read_csv("all12_events.csv")
    TABLES.mkdir(parents=True, exist_ok=True)
    lines = [
        r"\begin{table*}",
        r"\centering",
        r"\caption{Event-by-event reassessment of the 12 stronger COR2 anchors\label{tab:events}}",
        r"\scriptsize",
        r"\setlength{\tabcolsep}{3.1pt}",
        r"\begin{tabular}{rllrrrrll}",
        r"\toprule",
        r"No. & UTC & $N_{\rm pk}$ & $p_{\rm cell}$ & $p_{\rm best}$ & $p_{\rm common}$ & $|\Delta t_{tB-pB}|$ & Kinematic family & Class \\",
        r"\midrule",
    ]
    for i, r in enumerate(rows):
        utc = r["utc"].replace("2008-", "").replace("T", " ").replace(":30Z", "")
        kin = (r["best_model_family"] + " "
               + f"{float(r['best_v_inner_km_s']):.0f}"
               + r"$\to$"
               + f"{float(r['best_v_outer_km_s']):.0f}")
        end = r" \\"
        lines.append(
            f"{r['event_number']} & {utc} & {r['strict_matched_peak_count']} & "
            f"{float(r['strict_cell_stat_p_raw']):.3f} & {float(r['best_ridge_p_raw']):.3f} & "
            f"{float(r['common_ridge_p_raw']):.3f} & {float(r['tb_pb_global_phase_difference_min']):.0f} min & "
            f"{kin} & {r['category']}{end}"
        )
    lines += [
        r"\bottomrule",
        r"\end{tabular}",
        r"\begin{minipage}{0.98\textwidth}\vspace{2pt}\footnotesize\textit{Note.} All probabilities are raw empirical values. The kinematic family is a frozen alignment-model preference, not a bulk-flow measurement. Class A denotes expansion-stable outward ordering; B, coronal-compatible multi-ridge morphology; C, weak multi-ridge morphology; and D, an isolated or incomplete chain.\end{minipage}",
        r"\end{table*}",
    ]
    (TABLES / "all12_events.tex").write_text("\n".join(lines) + "\n")


def write_event_physics_table() -> None:
    rows = [
        ("1", "4 matched peaks", "30 min", "B", "Best non-pilot common morphology; no source clock or jump closure"),
        ("2", "2 peaks", "75 min", "D", "Too short for a train"),
        ("3", "1 peak", "75 min", "D", "Isolated enhancement"),
        ("4", "3 nearly regular peaks", "60 min", "B", "Cell-like hint; product timing disagrees"),
        ("5", "1 peak", "90 min", "D", "A timing hint without a radial chain"),
        ("6", "3 weak peaks", "180 min", "C", "Products describe different timing"),
        ("7", "1 peak", "45 min", "D", "No chain despite a fitted accelerating path"),
        ("8", "3 irregular peaks", "15 min", "C", "Geometry compatible, but weak spatial order"),
        ("9", "4 ordered nodes", "30 min", "A", "Strongest expansion-stable spatial sequence"),
        ("10", "3 regular peaks", "135 min", "B", "Morphology without preserved event phase"),
        ("11", "3 peaks", "30 min", "B", "Compatible accelerating-pattern candidate"),
        ("12", "3 ordered nodes", "0 min", "A", "Clearest individual ridge; possible gate-height coincidence"),
    ]
    lines = [
        r"\begin{table*}", r"\centering",
        r"\caption{Physical reading of every event in the primary sample\label{tab:eventphysics}}",
        r"\scriptsize", r"\setlength{\tabcolsep}{4.2pt}",
        r"\begin{tabular}{clccl}", r"\toprule",
        r"No. & Coronal morphology & $tB$--$pB$ offset & Class & Event-level inference \\",
        r"\midrule",
    ]
    for n, morph, phase, cls, inference in rows:
        lines.append(f"{n} & {morph} & {phase} & {cls} & {inference}" + r" \\")
    lines += [
        r"\bottomrule", r"\end{tabular}",
        r"\begin{minipage}{0.98\textwidth}\vspace{2pt}\footnotesize\textit{Note.} Class A means expansion-stable outward order; B, coronal-compatible multi-ridge morphology; C, weak multi-ridge morphology; and D, an incomplete chain.  The labels describe what the images support; none identifies an MHD branch.\end{minipage}",
        r"\end{table*}",
    ]
    (TABLES / "all12_event_physics.tex").write_text("\n".join(lines) + "\n")


def main() -> None:
    style()
    figure1_concept()
    figure2_survey()
    figure3_maps()
    figure4_expansion()
    figure5_height_controls()
    figure6_synthesis()
    figure7_resonator_height_chain()
    write_event_table()
    write_event_physics_table()
    print("Wrote seven vector PDF figures and two LaTeX event tables.")


if __name__ == "__main__":
    main()
