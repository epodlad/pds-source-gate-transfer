"""Stage verified release inputs in a fresh tree and run unmodified COR1 analysis programs."""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
META = ROOT / "data/cor1-release"
PROGRAMS = [
    "analyze_habbal_standing_mhd_branch", "analyze_node_activation_order_20080113",
    "analyze_dynamic_geometry_gate", "analyze_event09_2d_xfront_20080113",
    "analyze_two_candidate_2d_xfront_20080113", "analyze_two_event_representation_robustness_20080113",
    "analyze_expanding_front_sensitivity_20080113", "analyze_spatially_expanding_ridge_20080113",
]


def compare(a, b, path=""):
    if isinstance(a, dict):
        if set(a) != set(b):
            raise ValueError(f"Different keys: {path}")
        for k in a:
            compare(a[k], b[k], path + "/" + k)
    elif isinstance(a, list):
        if len(a) != len(b):
            raise ValueError(f"Different lengths: {path}")
        for i,(x,y) in enumerate(zip(a,b)):
            compare(x,y,f"{path}/{i}")
    elif isinstance(a,float) and isinstance(b,(int,float)):
        if not (math.isclose(a,b,abs_tol=1e-10,rel_tol=1e-8) or math.isnan(a) and math.isnan(b)):
            raise ValueError(f"{path}: {a} vs {b}")
    elif a != b:
        raise ValueError(f"{path}: {a!r} vs {b!r}")


def run(assets, work, all_tests):
    work.mkdir(parents=True,exist_ok=False)
    with (META / 'manifest.csv').open() as fp:
        entries={row['name']: row for row in csv.DictReader(fp)}
    mapping=json.loads((META/'input_paths.json').read_text())
    for name, rel in mapping.items():
        source=assets/name
        with source.open('rb') as stream:
            if hashlib.file_digest(stream,'sha256').hexdigest()!=entries[name]['sha256']:
                raise ValueError(f"Input checksum mismatch: {name}")
        target=work/rel
        target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy(source,target)
    for p in (ROOT/'analysis/frozen_programs').glob('*.py'):
        shutil.copy(p,work/p.name)
    env=dict(os.environ,MPLBACKEND='Agg',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
    for name in PROGRAMS if all_tests else PROGRAMS[:2]:
        print('Running',name,flush=True)
        with (work/(name+'.log')).open('w')as log:
            subprocess.run([sys.executable,str(work/(name+'.py'))],cwd=work,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
    pairs=[('habbal_standing_mhd_branch_results_no_BH.json',ROOT/'data/habbal_standing_mhd_results.json'),
           ('pds_20080113_node_activation_order_results_no_BH.json',assets/'pds_20080113_node_activation_order_results_no_BH.json')]
    if all_tests:
        for n in ('dynamic_geometry_gate_results_no_BH.json','pds_20080113_event09_2d_xfront_results_no_BH.json','pds_20080113_two_candidate_2d_xfront_results_no_BH.json','pds_20080113_two_event_representation_robustness_results.json','pds_20080113_expanding_front_sensitivity_results_no_BH.json','pds_20080113_spatially_expanding_ridge_results_no_BH.json'):
            pairs.append((n,assets/n))
    for name,ref in pairs:
        matches=list((work/'results').rglob(name))
        if len(matches)!=1:
            raise ValueError(f"Expected one output named {name}, found {len(matches)}")
        compare(json.loads(ref.read_text()),json.loads(matches[0].read_text()))
        print('MATCH',name,flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('assets',type=Path)
    p.add_argument('work',type=Path,help='New directory for this run')
    p.add_argument('--all',action='store_true')
    args=p.parse_args()
    run(args.assets.resolve(),args.work.resolve(),args.all)
