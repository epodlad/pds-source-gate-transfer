#!/usr/bin/env python3
"""Build a no-recalculation coronal reassessment of all 12 January events."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "results" / "pds_20080111_14_nonlinear_level0"
OUT = SOURCE / "all12_coronal_reassessment"
OUT.mkdir(parents=True, exist_ok=True)

PHASE_FILE = SOURCE / "event_phase_jitter" / \
    "pds_20080111_14_event_phase_jitter_results.json"
CELL_FILE = SOURCE / "event12_shock_cells" / \
    "pds_20080113_2137_event12_shock_cells_results.json"
GEOMETRY_FILE = ROOT / "results" / "dynamic_geometry_gate" / \
    "dynamic_geometry_gate_results_no_BH.json"
PILOT_FILE = SOURCE / "spatially_expanding_ridge_pilot" / \
    "pds_20080113_spatially_expanding_ridge_results_no_BH.json"
TIMING_FILE = SOURCE / "expanding_front_sensitivity" / \
    "pds_20080113_expanding_front_sensitivity_results_no_BH.json"
MATRIX_FILE = ROOT / "results" / "pds_combined" / "pds_event_height_matrix.json"

CATEGORY_LABELS = {
    "A": "expansion-stable ordered candidate",
    "B": "multi-ridge coronal-compatible morphology",
    "C": "weak multi-ridge morphology",
    "D": "isolated or incomplete chain",
}


def finite(value):
    return value is not None and isinstance(value, (int, float)) and math.isfinite(value)


def phase_label(difference):
    value = abs(float(difference))
    if value <= 30.0:
        return "timing-compatible within 0-30 min; source clock not established"
    if value <= 60.0:
        return "marginal after 15-30 min broadening; source clock not established"
    return "tB/pB timing unresolved after bounded broadening; source clock not established"


def classify(event_number, peaks, p_cell, p_best, p_common, pilot):
    if event_number in (9, 12):
        name = f"event{event_number}"
        ordered_reps = sum(
            pilot["events"]["measured"][name][rep]["p_outward_order"] >= 0.50
            for rep in ("bff", "base60", "nrgf60")
        )
        if ordered_reps >= 2 and pilot["pair_null"]["measured"][
            "primary_pair_p_raw"
        ] < 0.05:
            return "A"
    screening = min(p_cell, p_best, p_common)
    if peaks >= 3 and screening <= 0.10:
        return "B"
    if peaks >= 3:
        return "C"
    return "D"


def geometry_label(event_number, geometry_row, timing):
    if event_number in (9, 12):
        name = f"event{event_number}"
        ratios = [
            node["area_ratio"]
            for node in timing["events"][name]["base60"]["nodes"]
        ]
        prefix = f"measured pB width proxy A/Ainner up to {max(ratios):.2f}"
        if event_number == 12:
            return prefix + "; exploratory pB geometry/axis coincidence near 2.70-2.75 R_sun"
        return prefix + "; expanding-ridge geometry is sensitivity-stable"
    ratio = geometry_row.get("area_expansion_ratio")
    radius = geometry_row.get("gradient_peak_radius_rsun")
    if finite(ratio) and finite(radius):
        return (
            f"measured projected expansion proxy={ratio:.2f} near "
            f"{radius:.2f} R_sun; geometry is not causal evidence"
        )
    return "measured width proxy unavailable; spherical r^2 bracket only"


def kinematic_label(category, event_number, phase_row, pilot):
    model = phase_row["best_model"]
    family = model["family"]
    v0 = model["v_inner_km_s"]
    v1 = model["v_outer_km_s"]
    frozen = f"frozen {family} model {v0:.0f}->{v1:.0f} km/s"
    if category == "A":
        name = f"event{event_number}"
        values = [
            pilot["events"]["measured"][name][rep]["p_transport_25_300"]
            for rep in ("bff", "base60", "nrgf60")
        ]
        return (
            f"{frozen}; Ptransport={min(values):.3f}-{max(values):.3f}; "
            "outward order supported, transport not closed"
        )
    if category == "D":
        return f"{frozen}; chain transport not estimable from fewer than three matched peaks"
    return (
        f"{frozen}; edge speeds 20/25/300 are not hard failures, but no "
        "event-specific transport posterior was computed"
    )


def reassessment_text(event_number, category, peaks, spacing_cv, p_cell, p_best,
                      p_common, phase_difference):
    if category == "A" and event_number == 9:
        return (
            "Four matched peaks plus expansion-stable outward ordering.  This is the "
            "strongest spatial sequence, but its inter-node speeds remain poorly constrained."
        )
    if category == "A" and event_number == 12:
        return (
            "Three matched peaks and an individually identified ridge are preserved by "
            "the expanding filter.  It is closer to outward transport than #9, but remains "
            "below the full transport threshold."
        )
    if category == "B":
        if event_number == 1:
            return (
                "Four matched peaks and a low common-ridge raw p give the strongest "
                "non-pilot coronal-compatible morphology; no event clock or shock closure."
            )
        if event_number == 4:
            return (
                "Three nearly equally spaced peaks survive as a regular-cell hint; the "
                "60-min tB/pB mismatch makes transport/clock identification marginal."
            )
        if event_number == 10:
            return (
                "Three regularly spaced peaks support morphology, while the 135-min global "
                "tB/pB mismatch argues against a preserved event phase."
            )
        if event_number == 11:
            return (
                "Three peaks, a borderline best-ridge raw p, and 30-min product agreement "
                "make this a compatible but unconfirmed accelerating-pattern candidate."
            )
        return "Multi-ridge morphology is compatible after coronal broadening but not detected."
    if category == "C":
        if event_number == 6:
            return (
                "Three peaks are present, but their raw statistics are weak and the products "
                "differ by 180 min; retain only a weak morphology label."
            )
        if event_number == 8:
            return (
                "Projected expansion and 15-min product agreement are compatible, but radial "
                "spacing is irregular and the multi-peak statistic is weak."
            )
        return "Several peaks exist, but the frozen statistics do not support a distinct chain."
    if event_number == 2:
        return "Two matched peaks are insufficient for a periodic or shock-cell chain."
    if event_number == 3:
        return "One matched peak is an isolated enhancement, not a cross-height chain."
    if event_number == 5:
        return (
            "A common-ridge hint exists, but only one strict matched peak and 90-min product "
            "mismatch prevent a radial-chain interpretation."
        )
    if event_number == 7:
        return (
            "Only one matched peak is present; the accelerating fit is a model preference, "
            "not evidence for transported PDS."
        )
    return "The frozen measurements do not establish a cross-height PDS chain."


def build_rows():
    phase_data = json.loads(PHASE_FILE.read_text(encoding="utf-8"))
    cell_data = json.loads(CELL_FILE.read_text(encoding="utf-8"))
    geometry_data = json.loads(GEOMETRY_FILE.read_text(encoding="utf-8"))
    pilot = json.loads(PILOT_FILE.read_text(encoding="utf-8"))
    timing = json.loads(TIMING_FILE.read_text(encoding="utf-8"))
    matrix = json.loads(MATRIX_FILE.read_text(encoding="utf-8"))

    phase_rows = {row["event_number"]: row for row in phase_data["events"]}
    cell_rows = {row["event_number"]: row for row in cell_data["events"]}
    geometry_rows = {
        row["event_number"]: row
        for row in geometry_data["event_table_primary_total_b_200_km_s"]
    }
    rows = []
    for number in range(1, 13):
        phase_row = phase_rows[number]
        cell_row = cell_rows[number]
        geometry_row = geometry_rows[number]
        strict = cell_row["configs"]["strict"]
        peaks = int(strict["matched_peak_count"])
        spacing_cv = strict["spacing_cv"]
        p_cell = float(cell_row["strict_cell_stat_p_raw"])
        p_count = float(cell_row["strict_count_p_raw"])
        p_best = float(phase_row["best_ridge_p_raw"])
        p_common = float(phase_row["common_ridge_p_raw"])
        category = classify(number, peaks, p_cell, p_best, p_common, pilot)
        difference = float(phase_row["global_phase_product_difference_min"])
        row = {
            "event_number": number,
            "utc": phase_row["utc"],
            "cor2_z": float(phase_row["cor2_z"]),
            "strict_matched_peak_count": peaks,
            "strict_matched_center_radii_rsun": "|".join(
                f"{value:.3f}" for value in strict["matched_center_radii_rsun"]
            ),
            "strict_spacing_cv": spacing_cv,
            "strict_count_p_raw": p_count,
            "strict_cell_stat_p_raw": p_cell,
            "best_ridge_p_raw": p_best,
            "common_ridge_p_raw": p_common,
            "best_model_family": phase_row["best_model"]["family"],
            "best_v_inner_km_s": float(phase_row["best_model"]["v_inner_km_s"]),
            "best_v_outer_km_s": float(phase_row["best_model"]["v_outer_km_s"]),
            "tb_pb_global_phase_difference_min": difference,
            "phase_clock_reassessment": phase_label(difference),
            "geometry_reassessment": geometry_label(number, geometry_row, timing),
            "kinematic_reassessment": kinematic_label(
                category, number, phase_row, pilot
            ),
            "category": category,
            "category_label": CATEGORY_LABELS[category],
            "coronal_reassessment": reassessment_text(
                number, category, peaks, spacing_cv, p_cell, p_best,
                p_common, difference,
            ),
            "shock_specific_evidence": "not established",
            "new_p_value_computed": False,
        }
        rows.append(row)

    jan_case = next(case for case in matrix["cases"] if case["case"] == "2008-01-11–14")
    context = {
        "aggregate_cor1_raw_p": jan_case["cor1_ensemble_p"],
        "high_cor1_2p5_3p0_raw_p": jan_case["height_p"][-1],
        "high_cor1_verdict": jan_case["verdict"],
        "measured_width_pair_raw_p": pilot["pair_null"]["measured"][
            "primary_pair_p_raw"
        ],
        "measured_width_pair_exceedances": pilot["pair_null"]["measured"][
            "primary_pair_exceedances"
        ],
        "measured_width_pair_controls": pilot["control_pair_count"],
        "geometry_dependence_detected": geometry_data["geometry_dependence_detected"],
        "dynamic_gate_detected": geometry_data["dynamic_gate_detected"],
        "periodic_gate_detected": geometry_data["periodic_gate_detected"],
        "mhd_branch_identified": geometry_data["mhd_branch_identified"],
    }
    return rows, context


def write_report(rows, context):
    counts = {category: sum(row["category"] == category for row in rows)
              for category in CATEGORY_LABELS}
    lines = [
        "# Coronal-physics reassessment of all 12 events (no ridge recalculation)",
        "",
        "## Overall result",
        "",
        "The coronal correction changes interpretation, not measurements.  Exact phase "
        "preservation, constant area, constant amplitude, and hard 25/300 km/s boundaries "
        "are no longer required.  No ridge or raw p value was changed.",
        "",
        f"- Category A, expansion-stable ordered candidates: {counts['A']}/12 (#9, #12).",
        f"- Category B, multi-ridge coronal-compatible morphology: {counts['B']}/12 "
        "(#1, #4, #10, #11).",
        f"- Category C, weak multi-ridge morphology: {counts['C']}/12 (#6, #8).",
        f"- Category D, isolated/incomplete chains: {counts['D']}/12 (#2, #3, #5, #7).",
        "",
        f"At the ensemble level, the previously frozen high-COR1 aggregate remains "
        f"raw p={context['aggregate_cor1_raw_p']:.4f}, with localization at 2.5--3.0 "
        f"R_sun raw p={context['high_cor1_2p5_3p0_raw_p']:.4f}.  The measured-width "
        f"#9 -> #12 pair has {context['measured_width_pair_exceedances']}/"
        f"{context['measured_width_pair_controls']} shifted exceedances, raw "
        f"p={context['measured_width_pair_raw_p']:.4f}.",
        "",
        "These results support intermittent outward-ordered morphology in a subset of "
        "events, not a universal phase-coherent clock, universal stationary shock, or "
        "closed ballistic transport chain.",
        "",
        "## Event-by-event reassessment",
        "",
        "| # | UTC | Peaks | Cell raw p | Best/common ridge raw p | Model (km/s) | "
        "tB-pB phase | Tier | Coronal interpretation |",
        "| ---: | --- | ---: | ---: | ---: | --- | ---: | --- | --- |",
    ]
    for row in rows:
        lines.append(
            f"| {row['event_number']} | {row['utc']} | "
            f"{row['strict_matched_peak_count']} | {row['strict_cell_stat_p_raw']:.4f} | "
            f"{row['best_ridge_p_raw']:.4f}/{row['common_ridge_p_raw']:.4f} | "
            f"{row['best_model_family']} "
            f"{row['best_v_inner_km_s']:.0f}->{row['best_v_outer_km_s']:.0f} | "
            f"{row['tb_pb_global_phase_difference_min']:.0f} min | "
            f"{row['category']} | {row['coronal_reassessment']} |"
        )
    lines += [
        "",
        "## Meaning of the correction",
        "",
        "- A strict failure caused only by a boundary speed or 15--30 min timing offset is "
        "now treated as compatibility/inconclusive, not as a physical rejection.",
        "- A failure caused by fewer than three cross-height peaks is not rescued by "
        "expansion or turbulence.",
        "- Phase compatibility is not evidence for a resonator; the event-specific source "
        "clock remains unestablished in all 12 events.",
        "- Projected width and area are geometry diagnostics, not magnetic flux-tube area "
        "measurements.",
        "- Shock evidence remains unestablished for every event because no event has the "
        "required plasma state and jump-condition closure.",
        "",
        "## Global physical picture",
        "",
        "The safest synthesis is an intermittent, expanding or reforming cusp/current-sheet "
        "gate.  It can produce outward-ordered density/morphology sequences in selected "
        "events while losing exact phase and speed coherence through expansion, line-of-sight "
        "integration, and small-scale coronal variability.  The data do not support one rigid "
        "pulse, one universal 120-min phase, or one universal stationary Habbal shock.",
        "",
        "The dynamic-geometry, periodic-gate, and MHD-branch tests remain negative at the "
        "declared level.  Event #12 retains only an exploratory geometry-controlled-gate "
        "coincidence near 2.70--2.75 R_sun.",
        "",
        "All p values in this report are copied from frozen analyses and are raw; no BH "
        "correction and no new p-value calculation were applied.",
    ]
    path = OUT / "PDS_20080111_14_all12_coronal_reassessment_no_recalculation.md"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def main():
    rows, context = build_rows()
    csv_path = OUT / "pds_20080111_14_all12_coronal_reassessment.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    report_path = write_report(rows, context)
    json_path = OUT / "pds_20080111_14_all12_coronal_reassessment.json"
    payload = {
        "analysis_status": "post-extraction coronal-physics reassessment",
        "new_ridge_extraction": False,
        "new_p_values_computed": False,
        "no_BH": True,
        "categories": CATEGORY_LABELS,
        "context": context,
        "events": rows,
        "verdict": (
            "Two expansion-stable ordered candidates (#9, #12), four additional "
            "multi-ridge compatible morphologies (#1, #4, #10, #11), two weak "
            "multi-ridge cases (#6, #8), and four incomplete chains (#2, #3, #5, #7). "
            "No event-specific clock, closed transport chain, or MHD shock is established."
        ),
        "files": {
            "report": report_path.name,
            "table": csv_path.name,
        },
    }
    json_path.write_text(json.dumps(payload, indent=2, allow_nan=True), encoding="utf-8")
    print(json.dumps({
        "counts": {category: sum(row["category"] == category for row in rows)
                   for category in CATEGORY_LABELS},
        "events": {row["event_number"]: row["category"] for row in rows},
        "context": context,
        "files": [str(report_path), str(csv_path), str(json_path)],
    }, indent=2))


if __name__ == "__main__":
    main()
