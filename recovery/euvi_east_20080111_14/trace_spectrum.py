"""Recover the plotted path on its native code-angle axis and test the released cubes.
Usage: python trace_spectrum.py /path/to/cubes /path/to/output
Keep the two historical PNGs alongside this script. No fit to the reference period.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
from scipy.signal import detrend, periodogram, lombscargle
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("cubes", type=Path)
parser.add_argument("output", type=Path)
args = parser.parse_args()
args.output.mkdir(parents=True, exist_ok=True)

# Centers of the tick marks in the original 2160 x 2700 PNG (pixel origin at top left).
ytick = np.array([530.5, 450.5, 371.5, 292.5, 212.5, 133.5, 54.5])
r_ticks = np.arange(11,18) / 10
xtick = {171: [122.5,378.5,634.5,890.5,1146.5,1402.5,1658.5,1914.5],
         195: [138.5,387.5,635.5,884.5,1133.5,1382.5,1631.5,1880.5]}
bounds = (1.10, 1.20)
def csv_write(path, columns, rows):
    with path.open("w",newline="") as dest:
        writer=csv.writer(dest)
        writer.writerow(columns)
        writer.writerows(rows)

fig, ax = plt.subplots(2, 3, figsize=(13.3,7.1), layout="constrained")
summary = {}
for i, wave in enumerate((171,195)):
    png = HERE / f"euvi{wave}_lower_resonator_test.png"
    image=np.asarray(Image.open(png).convert("RGB"))
    assert image.shape[:2]==(2700,2160)
    cube_file=args.cubes / f"euvi{wave}_reconstructed_sector.npz"
    with np.load(cube_file) as data:
        rr=data["radii_rsun"]; angles=data["legacy_pa_deg"]; cube=data["intensity_dn_per_s"]
        times=data["time_utc"].astype("datetime64[ms]")
    band=(rr>=bounds[0]-1e-9)&(rr<=bounds[1]+1e-9)
    radii=rr[band]
    xp=np.polyfit(np.arange(100,241,20),xtick[wave],1)
    yp=np.polyfit(r_ticks,ytick,1)
    col_lo,col_hi=np.rint(np.polyval(xp,[160,195])).astype(int)
    center=[]; spans=[]
    for row in np.rint(np.polyval(yp,radii)).astype(int):
        white=np.flatnonzero(image[row,col_lo:col_hi].min(axis=1)>245)+col_lo
        assert len(white)>0 and np.all(np.diff(white)==1)
        center.append((np.median(white)-xp[1])/xp[0]);spans.append(np.ptp(white)/xp[0])
    center=np.array(center)
    csv_write(args.output/f"path_{wave}.csv",["radius_Rsun","code_angle_deg","projected_line_span_deg"],zip(radii,center,spans))
    stamp=(times-times[0])/np.timedelta64(1,"m")
    slot=np.rint(stamp/10).astype(int); assert len(np.unique(slot))==len(slot)
    regular_time=np.arange(slot[-1]+1)*10.
    variants={}; all_series={}
    for offset in (-0.25,0.,0.25):
        # Equal weight per radius; within each radius average pixels in +/-2.5 degrees.
        mask=np.abs(angles[None,:]-center[:,None]-offset)<=2.5
        profiles=np.column_stack([np.nanmean(cube[:,j,mask[k]],axis=1,dtype=float)
                                  for k,j in enumerate(np.flatnonzero(band))])
        value=np.mean(profiles,axis=1); assert np.isfinite(value).all()
        regular=np.interp(regular_time,slot*10.,value)
        all_series[offset]=value
        spectra={}
        for window in ("hann","boxcar"):
            freq,power=periodogram(regular,fs=0.1,detrend="linear",window=window)
            valid=(freq>=1/130)&(freq<=1/80)
            k=np.flatnonzero(valid)[power[valid].argmax()]
            spectra[window]={"in_band_peak_minutes":1/freq[k],"in_band_fraction_of_non_dc_power":power[valid].sum()/power[1:].sum()}
            if offset==0 and window=="hann":
                positive=(freq>=1/240)&(freq<=1/45)
                ax[i,2].plot(1/freq[positive],power[positive],color="#165978")
                ax[i,2].axvspan(80,130,color="#dce8eb",zorder=-1)
                ax[i,2].axvline(128.2143,ls=":",color="#8c502c")
        # Least-squares spectrum on actual UTC samples; same frequency bins, no gap interpolation.
        freq_band=freq[valid]
        y=value-np.polyval(np.polyfit(stamp,value,1),stamp)
        lp=lombscargle(stamp,y,2*np.pi*freq_band,normalize=True)
        spectra["uneven_time_in_band_peak_minutes"]=1/freq_band[np.argmax(lp)]
        variants[str(offset)]=spectra
    csv_write(args.output/f"brightness_{wave}.csv",["utc","mean_dn_per_s"],zip(times.astype(str),all_series[0.]))
    # Colorbar inversion on the archived static map. Spearman correlation is insensitive to monotone scaling.
    cbx={171:2040,195:2006}[wave]
    palette=image[55:600,cbx,:].astype(float)
    from scipy.spatial import cKDTree
    tree=cKDTree(palette)
    aidx=np.flatnonzero((angles>=143)&(angles<=197))[::4]
    typical=np.median(cube[:,band,:],axis=0)
    yy=np.rint(np.polyval(yp,radii)).astype(int)
    xx=np.rint(np.polyval(xp,angles[aidx])).astype(int)
    pixels=image[yy[:,None],xx[None,:]]
    dist,idx=tree.query(pixels)
    # Mask overplotted lines and pixels with colors inconsistent with the archived colorbar.
    good=(dist<4)&(np.abs(angles[aidx][None,:]-center[:,None])>0.75)&np.isfinite(typical[:,aidx])
    rho=spearmanr(-idx[good],typical[:,aidx][good]).statistic
    summary[str(wave)]={"frames":len(times),"missing_10min_slots":sorted(set(range(len(regular_time)))-set(slot.tolist())),
       "cube_sha256":hashlib.file_digest(cube_file.open("rb"),"sha256").hexdigest(),
       "path_deg_range":[float(center.min()),float(center.max())],"max_projected_line_span_deg":max(spans),
       "archived_image_vs_cube_median_rank_rho":rho,"image_comparison_samples":int(good.sum()),"variants":variants}
    ax[i,0].imshow(typical,extent=[angles[0]-.125,angles[-1]+.125,radii[0]-.0025,radii[-1]+.0025],origin="lower",aspect="auto",cmap="magma",vmin=np.nanpercentile(typical,5),vmax=np.nanpercentile(typical,95))
    ax[i,0].plot(center,radii,color="white",lw=1.4)
    ax[i,0].set(xlabel="Archived code angle (deg)",ylabel="Radius ($R_\\odot$)",title=f"EUVI {wave} Å · recovered path")
    ax[i,1].plot(stamp/60,all_series[0.],lw=.8,color="#165978")
    ax[i,1].set(xlabel="Hours since first exposure",ylabel="Mean intensity (DN/s)",title="Brightness along the path")
    ax[i,2].set(xlabel="Period (min)",ylabel="PSD ((DN/s)² min)",title="Linear detrend · Hann window")
fig.savefig(args.output/"euvi_path_spectra.png",dpi=160)
(args.output/"diagnostics.json").write_text(json.dumps(summary,indent=2)+"\n")
print(json.dumps(summary,indent=2))
