"""Test the digitized EUVI trace against the full-frame release.
python full_frame_check.py VERIFIED_FRAMES.csv FITS_ZIP_DIRECTORY OUTPUT
The code reads the exact mission bytes recorded by the release and rejects missing frames.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
from zipfile import ZipFile

import numpy as np
from scipy.ndimage import map_coordinates
from scipy.signal import periodogram, detrend, lfilter, windows
from scipy.stats import spearmanr

HERE=Path(__file__).resolve().parent
parser=argparse.ArgumentParser(description=__doc__)
for name in ("manifest","archives","output"):
    parser.add_argument(name,type=Path)
parser.add_argument("--save-cubes", action="store_true", help="Save east-limb diagnostic cubes for distribution")
args=parser.parse_args()
args.output.mkdir(parents=True,exist_ok=True)
manifest={r["filename"]:r for r in csv.DictReader(args.manifest.open())}
assert len(manifest)==716
seen=set()
def read_frame(b):
    hdr={}
    for off in range(0,len(b),80):
        s=b[off:off+80].decode("ascii")
        if s[:8].strip()=="END": break
        if s[8:10]=="= ": hdr[s[:8].strip()]=s[10:].split("/")[0].strip().strip("'").strip()
    assert int(hdr['BITPIX'])==16 and hdr['CTYPE1']=='HPLN-TAN' and hdr['CTYPE2']=='HPLT-TAN'
    assert hdr['CUNIT1']==hdr['CUNIT2']=='arcsec' and float(hdr['CRVAL1'])==float(hdr['CRVAL2'])==0
    arr=np.frombuffer(b,offset=((off+80+2879)//2880)*2880,dtype='>i2',count=int(hdr['NAXIS1'])*int(hdr['NAXIS2'])).reshape(int(hdr['NAXIS2']),-1)
    return hdr,arr

def sample(img,h,r,angle):
    scale=np.pi/648000.
    rho=r*float(h['RSUN'])*scale
    lon=rho*np.sin(np.deg2rad(angle)); lat=rho*np.cos(np.deg2rad(angle))
    u=np.tan(lon)/scale/float(h['CDELT1'])
    v=np.tan(lat)/np.cos(lon)/scale/float(h['CDELT2'])
    pc=np.array([[float(h[f'PC{i}_{j}']) for j in (1,2)] for i in (1,2)])
    xy=np.linalg.solve(pc,np.array([u.ravel(),v.ravel()]))
    x=xy[0]+float(h['CRPIX1'])-1; y=xy[1]+float(h['CRPIX2'])-1
    raw=map_coordinates(img,[y,x],order=1,output=np.float64,mode='constant',cval=np.nan)
    return ((raw*float(h['BSCALE'])+float(h['BZERO']))/float(h['EXPTIME'])).reshape(u.shape)

def power_fraction(y):
    arr=detrend(y,axis=-1,type="linear")*windows.hann(y.shape[-1],sym=False)
    p=np.abs(np.fft.rfft(arr,axis=-1))**2
    p[...,1:-1]*=2
    freq=np.fft.rfftfreq(y.shape[-1],d=10.)
    band=(freq>=1/130)&(freq<=1/80)
    return p[...,band].sum(-1)/p[...,1:].sum(-1),freq[band],p[...,band]

def analyze(times,brightness):
    t=np.array(times,dtype="datetime64[ms]")
    minutes=(t-t[0])/np.timedelta64(1,"m"); slots=np.rint(minutes/10).astype(int)
    assert len(slots)==358 and slots[-1]==359 and np.all(np.diff(slots)>0)
    grid=np.arange(360)
    y=np.interp(grid,slots,brightness)
    score,f,p=power_fraction(y)
    adjacent=np.diff(slots)==1
    dev=brightness-np.polyval(np.polyfit(minutes,brightness,1),minutes)
    phi=np.corrcoef(dev[:-1][adjacent],dev[1:][adjacent])[0,1]
    assert abs(phi)<1
    rng=np.random.default_rng(20080111)
    exceed=0
    for _ in range(100):
        eps=rng.standard_normal((200,360));eps[:,0]/=np.sqrt(1-phi**2)
        noise=lfilter([1.],[1.,-phi],eps,axis=-1)
        for j in np.setdiff1d(grid,slots):noise[:,j]=(noise[:,j-1]+noise[:,j+1])/2
        null,_,_=power_fraction(noise)
        exceed+=int((null>=score).sum())
    return dict(in_band_peak_min=1/f[np.argmax(p)],power_fraction=float(score),fitted_phi=phi,
                n_sim=20000,exceedances=exceed,conditional_ar1_p=(1+exceed)/20001,
                scope="New statistic and fitted stationary Gaussian AR(1) null; conditional on ROI and null parameters.")

summary={}
for wave in (171,195):
    path=np.genfromtxt(HERE/f"path_{wave}.csv",delimiter=',',names=True)
    radii=path['radius_Rsun']; center=path['code_angle_deg']
    rr=radii[:,None]; offsets=np.arange(-2.5,2.501,.25)
    angles=90.-(center[:,None]+offsets)
    cal=np.load(HERE/f"archived_image_{wave}.npz")
    angles_check=cal['angles']; rr_check=cal['radii'][:,None]
    transforms=[(sign,shift) for sign in (1,-1) for shift in (0,90,180,270)]
    series=[]; static=[]; samples=[]; sample_names=[]
    radial_grid=np.linspace(1.0,1.70,71)
    angle_grid=np.linspace(100.,240.,281)
    zips=sorted(args.archives.glob(f"stereo_a_euvi_{wave}_*_fits.zip"))
    assert len(zips)==3
    for zfile in zips:
        with ZipFile(zfile) as arc:
            for member in sorted(arc.namelist()):
                name=Path(member).name
                if not name.endswith(('.fts','.fits')): continue
                assert name not in seen
                b=arc.read(member); ref=manifest[name]
                assert hashlib.sha256(b).hexdigest()==ref['sha256']
                hdr,img=read_frame(b)
                assert int(hdr['WAVELNTH'])==wave==int(ref['channel']) and hdr['DATE-OBS']==ref['date_obs']
                obs=sample(img,hdr,rr,angles)
                assert np.isfinite(obs).all()
                if args.save_cubes:
                    samples.append(sample(img,hdr,radial_grid[:,None],90.-angle_grid).astype("float32"))
                    sample_names.append(name)
                series.append((hdr['DATE-OBS'],name,obs.mean()))
                static.append(np.array([sample(img,hdr,rr_check,sgn*angles_check+shift)
                                        for sgn,shift in transforms],dtype='float32'))
                seen.add(name)
        print("Processed",zfile.name,flush=True)
    series.sort()
    if args.save_cubes:
        lookup={name:j for j,name in enumerate(sample_names)}
        order=[lookup[s[1]] for s in series]
        np.savez_compressed(args.output/f"euvi{wave}_east_limb_diagnostic.npz",
                            intensity_dn_per_s=np.asarray(samples)[order],
                            radius_rsun=radial_grid,archived_code_angle_deg=angle_grid,
                            helioprojective_angle_deg=90.-angle_grid,
                            time_utc=np.array([s[0] for s in series]),filename=np.array([s[1] for s in series]))
    with (args.output/f"east_brightness_{wave}.csv").open("w",newline="") as dest:
        w=csv.writer(dest);w.writerow(["utc","filename","mean_dn_per_s"]);w.writerows(series)
    hist=cal['ranks']; keep=cal['valid']
    med=np.nanmedian(static,axis=0)
    spatial={}
    for (sgn,shift),mesh in zip(transforms,med):
        correlations=[]
        for hs,ok,curve in zip(hist,keep,mesh):
            mask=ok&np.isfinite(curve)
            if mask.sum()>65: correlations.append(spearmanr(hs[mask],curve[mask]).statistic)
        spatial[f"{sgn}*angle+{shift}"]=float(np.median(correlations))
    stats=analyze([s[0] for s in series],np.array([s[2] for s in series]))
    summary[str(wave)]=dict(frames=len(series),geometry=spatial,diagnostic=stats)
assert set(manifest)==seen
(args.output/'full_frame_results.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
