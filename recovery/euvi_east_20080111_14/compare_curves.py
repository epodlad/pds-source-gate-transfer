"""Compare time-domain shapes with the digitized historical middle panel, at zero lag."""
from pathlib import Path
import json
import numpy as np
from PIL import Image
from scipy.signal import butter, detrend, periodogram, sosfiltfilt
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path(__file__).resolve().parent
fig,axes=plt.subplots(2,2,figsize=(12,6.7),layout='constrained')
res={}
for i,w in enumerate([171,195]):
    img=np.asarray(Image.open(root/f'euvi{w}_lower_resonator_test.png').convert('RGB'))
    # Matplotlib tab:blue in the middle panel; exclude spines and tick labels.
    ys,xs=np.where(np.all(img[738:1288,:,:]==[31,119,180],axis=-1));ys+=738
    unique=np.unique(xs)
    px=(141.5,1785.5) if w==171 else (156.5,1754.5)
    hours=10+60*(unique-px[0])/(px[1]-px[0])
    yt=np.array([np.median(ys[xs==c]) for c in unique])
    if w==171: old=(1004.5-yt)*.2/(1095.5-1004.5)
    else: old=(1020.5-yt)*.05/(1104.5-1020.5)
    np.savetxt(root/f'historical_filtered_{w}.csv',np.c_[hours,old],delimiter=',',header='hours_from_2008_01_11_00UTC,plotted_component',comments='')
    data=np.genfromtxt(root/f'east_brightness_{w}.csv',delimiter=',',names=True,dtype=None,encoding='utf8')
    t=(np.array(data['utc'],dtype='datetime64[ms]')-np.datetime64('2008-01-11'))/np.timedelta64(1,'h')
    grid=t[0]+np.arange(360)/6
    light=np.interp(grid,t,data['mean_dn_per_s'])
    old_grid=np.interp(grid,hours,old,left=np.nan,right=np.nan)
    filt=sosfiltfilt(butter(3,[1/140,1/70],fs=1/10,btype='bandpass',output='sos'),light)
    mask=(grid>grid[0]+3)&(grid<grid[-1]-3)&np.isfinite(old_grid)
    def normalize(series): return (series-np.mean(series[mask]))/np.std(series[mask])
    r=np.corrcoef(old_grid[mask],filt[mask])[0,1]
    res[str(w)]={'zero_lag_correlation':r,'scope':'Fixed order-3 zero-phase Butterworth 70–140 min, edge exclusion 3 h, no lag fit. Original filter and panel height not verified.'}
    axes[i,0].plot(grid,normalize(old_grid),lw=1,color='#666',label='Archived plotted component')
    axes[i,0].plot(grid,normalize(filt),lw=1,color='#006b83',label='Re-extracted on east limb')
    axes[i,0].set(xlabel='Hours from 11 January 00:00 UTC',ylabel='Standard deviations',title=f'{w} Å · filtered time series')
    axes[i,0].legend(loc='upper left',fontsize=7)
    freq,power=periodogram(light,fs=.1,window='hann',detrend='linear')
    within=(freq>=1/200)&(freq<=1/60)
    axes[i,1].axvspan(80,130,color='#dbe9e9',label='80–130 min')
    axes[i,1].plot(1/freq[within],power[within],'.-',color='#006b83',lw=1)
    if w==171: axes[i,1].axvline(128.21,color='#a34f32',ls=':',label='Archived 171 Å peak')
    axes[i,1].set(xlabel='Period (min)',ylabel='PSD ((DN/s)² min)',title=f'{w} Å · eastern trace, Hann spectrum')
    axes[i,1].legend(fontsize=7)
fig.savefig(root/'east_limb_time_series.png',dpi=150)
(root/'curve_comparison.json').write_text(json.dumps(res,indent=2)+'\n')
print(res)
