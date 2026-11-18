"""New, explicitly labelled reconstruction; original selection/grid unverified.
Input: archive candidate inventory CSV. Retains full FITS; builds exposure-normalized
polar sectors using header WCS, without background subtraction or noise filtering.
"""
import csv,json,hashlib,time,urllib.request,concurrent.futures,argparse
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import map_coordinates

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--inventory',required=True);ap.add_argument('--output',required=True);ap.add_argument('--workers',type=int,default=6);a=ap.parse_args()
 out=Path(a.output);raw=out/'raw';raw.mkdir(parents=True,exist_ok=True);cache=out/'sectors';cache.mkdir(exist_ok=True)
 rows=list(csv.DictReader(open(a.inventory)))
 selected=[]
 for row in rows:
  n=row['filename'];tail=n[12:15]
  if '_n4euA.fts' in n and tail in ('030','100'):
   row['expected_channel']=195 if tail=='030' else 171;selected.append(row)
 # Explicit NEW grid; not recovered from original analysis.
 radii=np.linspace(1.0,1.5,101);pa=np.linspace(140,200,241)
 def run(row):
  n=row['filename'];p=raw/n;target=cache/(n+'.npz');meta=cache/(n+'.json')
  if target.exists() and meta.exists():return json.loads(meta.read_text())
  for attempt in range(3):
   try:
    if not p.exists():
     with urllib.request.urlopen(row['url'],timeout=60) as r:b=r.read()
     temp=p.with_suffix('.part');temp.write_bytes(b);temp.replace(p)
    with fits.open(p,memmap=False) as hd:
     hd.verify('exception');h=hd[0].header;image=hd[0].data.astype(float)
     wave=int(h['WAVELNTH']);assert wave==row['expected_channel'],f'Unexpected wavelength {wave}'
     exp=float(h['EXPTIME']);assert exp>0
     # Legacy PA convention: x=+r sin(PA), y=+r cos(PA).
     x=radii[:,None]*float(h['RSUN'])*np.sin(np.deg2rad(pa))[None,:]
     y=radii[:,None]*float(h['RSUN'])*np.cos(np.deg2rad(pa))[None,:]
     w=WCS(h,naxis=2);xx,yy=w.all_world2pix(x/3600,y/3600,0)
     vals=map_coordinates(image/exp,[yy,xx],order=1,mode='constant',cval=np.nan).astype('float32')
     result=dict(filename=n,url=row['url'],channel=wave,date_obs=h['DATE-OBS'],exptime=exp,shape=list(image.shape),sha256=hashlib.sha256(p.read_bytes()).hexdigest(),finite_fraction=float(np.isfinite(vals).mean()),status='ok')
    np.savez_compressed(target,intensity=vals);meta.write_text(json.dumps(result,indent=2));return result
   except Exception as e:
    if attempt==2:return dict(filename=n,url=row['url'],status='error',error=str(e))
    time.sleep(1)
 results=[]
 with concurrent.futures.ThreadPoolExecutor(max_workers=a.workers) as pool:
  for result in pool.map(run,selected):
   results.append(result)
   if len(results)%20==0:print('Processed',len(results),'/',len(selected),'errors',sum(r['status']!='ok' for r in results),flush=True)
   (out/'progress.json').write_text(json.dumps(results,indent=2))
 for wave in (171,195):
  rr=sorted([r for r in results if r['status']=='ok' and r['channel']==wave],key=lambda r:r['date_obs'])
  data=np.stack([np.load(cache/(r['filename']+'.npz'))['intensity'] for r in rr])
  np.savez_compressed(out/f'euvi{wave}_reconstructed_sector.npz',intensity_dn_per_s=data,radii_rsun=radii,legacy_pa_deg=pa,time_utc=np.array([r['date_obs'] for r in rr]),filenames=np.array([r['filename'] for r in rr]))
 print('COMPLETE',len(results),'errors',sum(r['status']!='ok' for r in results),flush=True)
if __name__=='__main__':main()
