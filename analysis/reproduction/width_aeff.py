"""Exploratory January COR1 width and isothermal effective-area proxy, 1.8--2.3 solar radii."""
import argparse
import csv
from pathlib import Path
import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.constants import Boltzmann, proton_mass

GM_SUN = 1.3271244e20   # m^3/s^2, IAU nominal solar mass parameter
R_SUN = 6.957e8        # m, IAU nominal solar radius


def width(offsets, profile, outer):
    y = gaussian_filter1d(profile, 4., mode="nearest")  # one degree on a 0.25 degree grid
    side = y[(abs(offsets) >= outer-4) & (abs(offsets) <= outer)]
    core = np.flatnonzero((abs(offsets) <= 3) & np.isfinite(y))
    if np.isfinite(side).sum() < 8 or not len(core):
        return np.nan
    base = np.nanmedian(side)
    noise = 1.4826 * np.nanmedian(abs(side-base))
    peak = core[np.argmax(y[core])]
    if y[peak]-base <= 2*noise:
        return np.nan
    half = (y[peak]+base)/2
    crossings = []
    for direction in (-1,1):
        i = peak
        while 0 <= i+direction < len(y) and y[i+direction] >= half:
            i += direction
        j = i+direction
        if not 0 <= j < len(y) or not np.isfinite(y[j]):
            return np.nan
        crossings.append(np.interp(half,[y[j],y[i]],[offsets[j],offsets[i]]))
    return crossings[1]-crossings[0]


def minimum(r, values):
    if not np.isfinite(values).any():
        return "", False
    j = np.nanargmin(values)
    bracket = 0 < j < len(r)-1 and np.isfinite(values[j-1]) and np.isfinite(values[j+1])
    return round(r[j],3), bool(bracket)


def run(assets, output, mu):
    output.mkdir(parents=True,exist_ok=True)
    ct2 = Boltzmann*1e6/(mu*proton_mass)
    with np.load(assets/'cor1a_secchi_background_fitpol_sector_15min.npz') as cube, np.load(assets/'pds_20080111_14_level0_nonlinear_transport_maps.npz') as maps:
        t,rad,pa = (cube[k] for k in ('minutes','radii','pa'))
        path=maps['traced_path_pa']
        assert len(rad)==len(path)
        bins=np.flatnonzero((rad>=1.8-1e-8)&(rad<=2.3+1e-8))
        r=rad[bins]
        summaries=[]; profiles=[]
        for label, use in [('60h',(t>=720)&(t<4320)),('all',np.ones(len(t),bool))]:
            for product in ('total_b','pb_signed'):
                image=np.nanmedian(cube[product][use].astype(float),axis=0)
                for outer in (12,16,20):
                    offsets=np.arange(-outer,outer+0.125,0.25)
                    angles=[]
                    for i in bins:
                        keep=np.isfinite(image[i])
                        if keep.sum()<3:
                            angles.append(np.nan)
                            continue
                        sampled=np.interp(path[i]+offsets,pa[keep],image[i,keep],left=np.nan,right=np.nan)
                        angles.append(width(offsets,sampled,outer))
                    w=r*np.deg2rad(angles)
                    # An axisymmetric circular cross section is a geometric assumption.
                    area=np.pi*w**2/4
                    aeff=area*np.exp(GM_SUN/(r*R_SUN*ct2))
                    count=int(np.isfinite(w).sum())
                    summaries.append([label,int(use.sum()),product,outer,count,len(r),count/len(r),*minimum(r,w),*minimum(r,aeff)])
                    for rr,angle,ww,aa,eff in zip(r,angles,w,area,aeff):
                        profiles.append([label,product,outer,round(rr,3),angle,ww,aa,eff])
    for name,header,rows in [
        ('width_summary','window frames product half_window_deg n_valid n_total valid_fraction sampled_width_min_rsun width_bracketed sampled_aeff_min_rsun aeff_bracketed',summaries),
        ('width_profiles','window product half_window_deg radius_rsun fwhm_deg width_rsun area_rsun2 aeff_rsun2',profiles)]:
        with (output/(name+'.csv')).open('w',newline='') as f:
            writer=csv.writer(f);writer.writerow(header.split());writer.writerows(rows)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('assets',type=Path)
    parser.add_argument('output',type=Path)
    parser.add_argument('--mu',type=float,default=0.6,help='Mean mass per particle in proton masses')
    args=parser.parse_args()
    run(args.assets,args.output,args.mu)
