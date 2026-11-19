import json, csv, pathlib
ROOT=pathlib.Path(__file__).resolve().parent

def bh(p):
 order=sorted(range(len(p)),key=lambda i:p[i]); out=[0.]*len(p); last=1.
 for rank in range(len(p),0,-1):
  i=order[rank-1];last=min(last,p[i]*len(p)/rank);out[i]=last
 return out
checks=[]
for path in sorted((ROOT/'reference').glob('*.json')):
 d=json.loads(path.read_text())
 families={'metrics':[('p_red_noise_raw','q_BH_red_noise_raw'),('p_red_noise_filter_matched','q_BH_red_noise_filter_matched'),('p_phase_circular_shift','q_BH_phase_circular_shift')], 'svd_modes':[('p_red_noise','q_red_noise_BH'),('p_phase_circular_shift','q_phase_BH')]}
 for group,pairs in families.items():
  rows=d[group]
  for pk,qk in pairs:
   computed=bh([r[pk] for r in rows]); error=max(abs(a-r[qk]) for a,r in zip(computed,rows))
   checks.append(dict(file=path.name,group=group,p_field=pk,n=len(rows),max_error=error,passed=error<1e-12))
assert all(x['passed'] for x in checks)
d=json.loads((ROOT/'reference/euvi171_lower_resonator_results.json').read_text())
r=next(x for x in d['metrics'] if x['metric']=='brightness_1.10_1.20R')
assert round(r['p_red_noise_raw'],4)==0.0349
assert round(r['p_red_noise_filter_matched'],4)==0.0279
assert round(r['dominant_period_raw_min'],2)==128.21
(ROOT/'audit.json').write_text(json.dumps({'scope':'Saved results only; no FITS, spectra or null trials rerun','BH_checks':checks,'preprint_page7_values':r},indent=2)+'\n')
print('PASS: 10 BH families (120 q-values); 3 rounded page-7 numbers match saved EUVI171 row.')
print('Maximum q error:',max(x['max_error'] for x in checks))
