import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
from pathlib import Path

from scipy.interpolate import UnivariateSpline
from scipy.optimize import curve_fit

#%%
def sens_pydis(obj_wave, obj_flux, std_wave, std_flux, std_wth):
    # Automatically exclude these obnoxious lines...
    balmer = np.array([6563, 4861, 4341], dtype='float')

    # down-sample (ds) the observed flux to the standard's bins
    obj_flux_ds = []
    obj_wave_ds = []
    std_flux_ds = []
    for i in range(len(std_wave)):
        rng = np.where((obj_wave >= std_wave[i] - std_wth[i] / 2.0) &
                       (obj_wave < std_wave[i] + std_wth[i] / 2.0))
        IsH = np.where((balmer >= std_wave[i] - std_wth[i] / 2.0) &
                       (balmer < std_wave[i] + std_wth[i] / 2.0))

        # does this bin contain observed spectra, and no Balmer line?
        if (len(rng[0]) > 1) and (len(IsH[0]) == 0):
            # obj_flux_ds.append(np.sum(obj_flux[rng]) / std_wth[i])
            obj_flux_ds.append( np.nanmean(obj_flux[rng]) )
            obj_wave_ds.append(std_wave[i])
            std_flux_ds.append(std_flux[i])

    # the ratio between the standard star flux and observed flux
    # has units like erg / counts
    ratio = np.abs(np.array(std_flux_ds, dtype='float') /
                   np.array(obj_flux_ds, dtype='float'))
    spl = UnivariateSpline(obj_wave_ds, ratio, ext=0, k=2 ,s=0.0025)
    sensfunc2 = spl(obj_wave)
    
#    plt.figure()
#    plt.plot(std_wave, std_flux, 'r', alpha=0.5, label='standard flux')
#    plt.xlabel('Wavelength')
#    plt.ylabel('Standard Star Flux')
#    plt.legend()
#    plt.show()

#    plt.figure()
#    plt.plot(obj_wave, obj_flux, 'k', label='observed counts')
#    plt.plot(obj_wave_ds, obj_flux_ds, 'bo',
#            label='downsample observed')
#    plt.xlabel('Wavelength')
#    plt.ylabel('Observed Counts/S')
#    plt.legend()
#    plt.show()

#    plt.figure()
#    plt.plot(obj_wave_ds, ratio, 'ko', label='sensfunc')
#    plt.plot(obj_wave, sensfunc2, label='interpolated sensfunc')
#    plt.xlabel('Wavelength')
#    plt.ylabel('log Sensfunc')
#    plt.legend()
#    plt.show()
#
#    plt.figure()
#    plt.plot(obj_wave, obj_flux*(10**sensfunc2),'k',
#             label='applied sensfunc')
#    plt.plot(std_wave, std_flux, 'ro', alpha=0.5, label='standard flux')
#    plt.xlabel('Wavelength')
#    plt.ylabel('Standard Star Flux')
#    plt.legend()
#    plt.show()
#    
    return sensfunc2
    

#%%
top = Path("/home/astro_02/AO2019-2/")
modelpath = Path(top/"uka2i.dat")
modelspec = pd.read_csv(modelpath, skiprows=1, delimiter=',',
                      names=["Wavelength", "flux", "?1", "?2", "?3", "?4"])
usp3 = UnivariateSpline(modelspec["Wavelength"], modelspec["flux"], k=3, s=0)
usp1 = UnivariateSpline(modelspec["Wavelength"], modelspec["flux"], k=1, s=0)

specpath = top/"2019-10-24/Specdata"
hd207673specpaths = list(specpath.glob("HD207673*"))
WISEAspecpaths = list(specpath.glob("WISEA_raw*.csv"))




#%%
senss1 = []
senss3 = []

fig, axs = plt.subplots(1,1)

for fpath in hd207673specpaths:
    stdobsspec = pd.read_csv(fpath)
#    axs.plot(stdobsspec["Wavelength"], stdobsspec["Value"])
    std1 = usp1(stdobsspec["Wavelength"])
    std3 = usp3(stdobsspec["Wavelength"])
    senss1.append(std1/stdobsspec["Value"])
    senss3.append(std3/stdobsspec["Value"])
    
    
senss1 = senss1[2]
senss3 = senss3[2]
#senss1 = np.median(senss1, axis=0)
#senss3 = np.median(senss3, axis=0)

#ax2.plot(stdobsspec["Wavelength"], senss1, 'k-')
#ax2.plot(stdobsspec["Wavelength"], senss3, 'g-')

#for fpath in WISEAspecpaths:
#    rawspec = pd.read_csv(fpath)
#    axs.plot(rawspec["Wavelength"], rawspec["Value"]*senss1, 'k-', lw=0.5)
#    axs.plot(rawspec["Wavelength"], rawspec["Value"]*senss3, 'g-', lw=0.5)
    

specs1 = []
for fpath in WISEAspecpaths:
    rawspec = pd.read_csv(fpath)
    spec1 = rawspec["Value"]*senss1
    spec3 = rawspec["Value"]*senss3
    specs1.append(spec1)
#    axs.plot(rawspec["Wavelength"],
#             bn.move_mean(rawspec["Value"]*senss1, window=5, min_count=1), 'k-')
#    axs.plot(rawspec["Wavelength"], spec1 , 'k-', lw=0.5)
#    axs.plot(rawspec["Wavelength"], spec3, 'g-', lw=0.5)
    
med = np.median(specs1, axis=0)
std = np.std(specs1, axis=0, ddof=1)
axs.set_xlabel(r'$\lambda$ ($\AA$)')
axs.set_ylabel('Arbitrary Flux')
axs.set_title("WISEA Spectrum")
axs.plot(rawspec["Wavelength"], med, 'k-', lw=0.5)
axs.plot(rawspec["Wavelength"], np.median(specs1, axis=0) , 'k-', lw=0.5)
axs.fill_between(rawspec["Wavelength"], med-std, med+std, alpha=0.5, color='r', lw=0)
    
    
ax2 = axs.twinx()
ax2.plot(rawspec["Wavelength"], senss1
         ,color='g', lw=1.5,alpha = 0.5)
#%%
sens2 = sens_pydis(stdobsspec["Wavelength"].to_numpy(), 
               stdobsspec["Value"].to_numpy(), 
               modelspec["Wavelength"].to_numpy(), 
               modelspec["flux"].to_numpy(), 
               std_wth=np.ones_like(modelspec["Wavelength"])*5)
#%%

fig, axs = plt.subplots(1,1)
specs=[]

for fpath in WISEAspecpaths:
    rawspec = pd.read_csv(fpath)
    specs.append(rawspec["Value"]*sens2)
    axs.plot(rawspec["Wavelength"], rawspec["Value"]*sens2, 'k-', lw=0.5)

#%%
specs
#%%
mspec=[]

for i in range(len(specs[0])):
    m = []
    for j in range(len(specs)):
        m.append(specs[j][i])
    mspec.append(np.median(m))
#%%
fig, axs = plt.subplots(1,1)


axs.plot(rawspec["Wavelength"], mspec, 'k-', lw=0.5)

