#!/usr/bin/env python3
"""Lightweight integrity checks for the compact PDS software archive."""
from __future__ import annotations

import ast
import csv
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
CODE = ROOT / "analysis" / "frozen_programs"

REQUIRED = [
    DATA / "all12_events.csv",
    DATA / "event_phase_jitter_table.csv",
    DATA / "event_phase_jitter_results.json",
    DATA / "height_matrix_no_BH.csv",
    DATA / "resonator_candidate_metrics.csv",
    DATA / "spatial_expanding_posteriors.npz",
    DATA / "two_event_maps.npz",
]


def main() -> None:
    missing = [str(p.relative_to(ROOT)) for p in REQUIRED if not p.exists()]
    if missing:
        raise SystemExit(f"Missing required files: {missing}")

    py_files = sorted(CODE.glob("*.py"))
    if not py_files:
        raise SystemExit("No frozen Python programs found")
    for path in py_files:
        ast.parse(path.read_text(encoding="utf-8"), filename=str(path))

    with (DATA / "all12_events.csv").open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 12, f"Expected 12 frozen events, found {len(rows)}"

    with (DATA / "event_phase_jitter_results.json").open(encoding="utf-8") as handle:
        phase = json.load(handle)
    assert phase.get("event_count") == 12
    assert phase.get("no_BH") is True

    with (DATA / "resonator_candidate_metrics.csv").open(newline="", encoding="utf-8") as handle:
        metrics = list(csv.DictReader(handle))
    assert any(row.get("domain") == "EUVI 171" for row in metrics)

    for name in ("spatial_expanding_posteriors.npz", "two_event_maps.npz"):
        with np.load(DATA / name, allow_pickle=False) as archive:
            assert archive.files, f"{name} contains no arrays"

    print(f"OK: {len(py_files)} frozen Python programs parsed successfully")
    print("OK: 12-event table and compact JSON/CSV/NPZ products validated")


if __name__ == "__main__":
    main()
