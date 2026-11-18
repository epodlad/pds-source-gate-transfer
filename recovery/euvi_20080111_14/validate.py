"""Validate rebuilt cube provenance and sampled exposure normalization."""
import json,sys,hashlib
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import map_coordinates
root=Path(sys.argv[1]);rows=json.loads((root/'progress.json').read_text());report={}
for wave in (171,195):
 p=root/f'euvi{wave}_reconstructed_sector.npz'
 with np.load(p,allow_pickle=False) as z:
  v=z['intensity_dn_per_s'];t=z['time_utc'];names=z['filenames'];r=z['radii_rsun'];pa=z['legacy_pa_deg']
  assert v.shape==(len(t),len(r),len(pa));assert len(set(names))==len(names)
  assert np.all(np.diff(t.astype('datetime64[ms]').astype('int64'))>0)
  rec={q['filename']:q for q in rows if q['status']=='ok' and q['channel']==wave};assert set(names)==set(rec)
  errors=[]
  for idx in sorted(set([0,len(t)//2,len(t)-1])):
   raw=root/'raw'/names[idx];assert hashlib.sha256(raw.read_bytes()).hexdigest()==rec[names[idx]]['sha256']
   with fits.open(raw,memmap=False) as hd:
    hd.verify('exception');h=hd[0].header;assert int(h['WAVELNTH'])==wave
    ri=np.array([10,40,80]);pi=np.array([40,120,200]);xx=r[ri]*float(h['RSUN'])*np.sin(np.deg2rad(pa[pi]));yy=r[ri]*float(h['RSUN'])*np.cos(np.deg2rad(pa[pi]))
    x,y=WCS(h,naxis=2).all_world2pix(xx/3600,yy/3600,0)
    expected=map_coordinates(hd[0].data.astype(float)/float(h['EXPTIME']),[y,x],order=1,mode='constant',cval=np.nan)
    got=v[idx,ri,pi];assert np.allclose(expected,got,rtol=1e-6,equal_nan=True)
    errors.append(float(np.nanmax(abs(expected-got))))
  report[str(wave)]={'frames':len(t),'shape':list(v.shape),'first':str(t[0]),'last':str(t[-1]),'finite_fraction':float(np.isfinite(v).mean()),'max_sample_abs_error':max(errors),'cube_sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
(root/'validation.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
