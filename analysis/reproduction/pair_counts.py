"""Recompute panel-b counts separately from the opposite-diagonal contrast test."""
import csv
import sys
from pathlib import Path

def rows(path):
    with path.open(newline='') as f:
        return list(csv.DictReader(f))

def flag(value):
    return {'True': True, 'False': False}[value]

def counts(controls, name, columns, archived):
    selected = []
    for row in controls:
        passed = sum(flag(row[c]) for c in columns) >= 2
        if passed != flag(row[archived]):
            raise ValueError(f"Pair flags disagree for {row['early_anchor_min']}")
        if passed:
            selected.append((float(row['early_anchor_min']), float(row['late_anchor_min'])))
    k,n=len(selected),len(controls)
    print(f'{name},{k},{n},{(k+1)/(n+1):.10f},{selected}')

def main(spatial, alternative=None):
    a=rows(spatial)
    keys=lambda r: (float(r['early_anchor_min']),float(r['late_anchor_min']))
    times=set(map(keys,a))
    assert len(times)==len(a) and all(l-e==360 for e,l in times)
    for law in ('fixed','spherical','measured'):
        cols=[f'{law}_{rep}_pair_exceed' for rep in ('bff','base60','nrgf60')]
        counts(a,law,cols,f'{law}_primary_pair_exceed')
    if alternative:
        b=rows(alternative)
        assert len(b)==len(a) and set(map(keys,b))==times
        assert all(flag(r['fully_valid']) for r in b)
        cols=[f'{rep}_pair_pass' for rep in ('base60','base120','nrgf60')]
        counts(b,'X-contrast',cols,'primary_pair_pass')
if __name__=='__main__':
    main(*(Path(a) for a in sys.argv[1:]))
