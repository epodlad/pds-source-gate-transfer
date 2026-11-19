"""Re-run the frozen phase analysis and figure generator from archived numerical inputs."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
MAP = "pds_20080111_14_level0_nonlinear_transport_maps.npz"
STEM = "pds_20080111_14_event_phase_jitter"


def run(script):
    subprocess.run([sys.executable, str(script)], check=True, env={**os.environ, "MPLBACKEND": "Agg"})


def reproduce(out):
    out.mkdir(parents=True, exist_ok=False)
    phase_dir = out / "results" / "pds_20080111_14_nonlinear_level0"
    phase_dir.mkdir(parents=True)
    shutil.copy2(HERE / MAP, phase_dir)
    shutil.copy2(REPO / "data/level0_nonlinear_transport_results.json",
                 phase_dir / "pds_20080111_14_level0_nonlinear_transport_results.json")
    for script in ("analyze_pds_event_phase_jitter_20080111_14.py",
                   "analyze_pds_nonlinear_transport_20080111_14.py"):
        shutil.copy2(REPO / "analysis/frozen_programs" / script, out)
    run(out / "analyze_pds_event_phase_jitter_20080111_14.py")
    result = json.loads((phase_dir / "event_phase_jitter" / (STEM + "_results.json")).read_text())
    reference = json.loads((REPO / "data/event_phase_jitter_results.json").read_text())
    if result != reference:
        raise ValueError("The recalculated event-phase JSON differs from the archived reference.")
    for name in ("analysis", "figures", "tables"):
        (out / name).mkdir()
    shutil.copytree(REPO / "data", out / "data")
    script = out / "analysis/make_manuscript_figures.py"
    shutil.copy2(HERE / script.name, script)
    run(script)
    print("Matched the archived phase results for all 12 events. Outputs:", out)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path, help="A new directory for the outputs")
    reproduce(parser.parse_args().output.resolve())
