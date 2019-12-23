import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
from pathlib import Path

from scipy.interpolate import UnivariateSpline
from scipy.optimize import curve_fit

#%%
def gauss1(x, a, mu, s):
    return a*np.exp(-(x-mu)**2/(2*s**2))


def gauss2(x, a1, a2, mu1, mu2, s1, s2, slope, const):
    return gauss1(x, a1, mu1, s1) + gauss1(x, a2, mu2, s2) + slope*x + const

def gaus32(x, a1, a2, a3, mu1, mu2, mu3, s1, s2, s3, slope, const):
    return gauss1(x, a1, mu1, s1) + gauss1(x, a2, mu2, s2) + gauss1(x, a3, mu3, s3)+ slope*x + const

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
#
#    plt.figure()
#    plt.plot(obj_wave, obj_flux, 'k', label='observed counts')
#    plt.plot(obj_wave_ds, obj_flux_ds, 'bo',
#            label='downsample observed')
#    plt.xlabel('Wavelength')
#    plt.ylabel('Observed Counts/S')
#    plt.legend()
#    plt.show()
#
#    plt.figure()
#    plt.plot(obj_wave_ds, ratio, 'ko', label='sensfunc')
#    plt.plot(obj_wave, sensfunc2, label='interpolated sensfunc')
#    plt.xlabel('Wavelength')
#    plt.ylabel('log Sensfunc')
#    plt.legend()
#    plt.show()
##
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
top = Path("/home/astro_02/AO2019-2")
modelpath = top/"uka0iii.dat"
modelspec = pd.read_csv(modelpath, skiprows=1, delimiter=',',
                        names=["Wavelength", "flux", "?1", "?2", "?3"])
usp3 = UnivariateSpline(modelspec["Wavelength"], modelspec["flux"], k=3, s=0)
usp1 = UnivariateSpline(modelspec["Wavelength"], modelspec["flux"], k=1, s=0)

specpath = top/"2019-10-24/Specdata"
stdspecpaths = list(specpath.glob("HD61931*"))
ngcspecpaths = list(specpath.glob("NGC2639_raw*.csv"))

#%%
senss1 = []
senss3 = []

for fpath in stdspecpaths:
    stdobsspec = pd.read_csv(fpath)
#    axs.plot(stdobsspec["Wavelength"], stdobsspec["Value"])
    std1 = usp1(stdobsspec["Wavelength"])
    std3 = usp3(stdobsspec["Wavelength"])
    senss1.append(std1/stdobsspec["Value"])
    senss3.append(std3/stdobsspec["Value"])

#%%
xtest=np.arange(len(senss1[0]))

ftsens=[]
ftzero_sens=[]
for i in range(len(senss1)):
    ftzero_sens.append(senss1[i][0])
    ftsens.append(senss1[i][1100]-senss1[i][0])
print( np.where(ftzero_sens==np.median(ftzero_sens)), np.where(ftsens==np.median(ftsens)))

#%%
fact = []
for i in range(len(senss1)):
    fact.append( ftsens[8]/ftsens[i] )
    senss1[i] = fact[i]*( senss1[i] - senss1[i][0] )+senss1[1][0]     

master_sens = []
for i in range(len(senss1[0])):
    x = []
    for j in range(len(senss1)):
        x.append(senss1[j][i])
    master_sens.append(np.median(x))
#%%
fig, axs = plt.subplots(1,1)
for i in range(len(senss1)):

    axs.plot(xtest, senss1[i],
#         color='k',
         label = i,
         lw=1.5,
         alpha = 0.2)
    
axs.plot(xtest, master_sens,
         color='k',
         lw=1.5,
         alpha = 1)

    
axs.legend()

#%%




fig, axs = plt.subplots(1,1)




senss1 = senss1[2]
senss3 = senss3[2]
#senss1 = np.median(senss1, axis=0)
#senss3 = np.median(senss3, axis=0)

ax2 = axs.twinx()
#ax2.plot(stdobsspec["Wavelength"], senss1, 'k-')
#ax2.plot(stdobsspec["Wavelength"], senss3, 'g-')

specs1 = []
for fpath in ngcspecpaths:
    rawspec = pd.read_csv(fpath)
    mspec = rawspec["Value"]*master_sens
    specs1.append(mspec)
#    spec1 = rawspec["Value"]*senss1
#    spec3 = rawspec["Value"]*senss3
#    specs1.append(spec1)
#    axs.plot(rawspec["Wavelength"],
#             bn.move_mean(rawspec["Value"]*senss1, window=5, min_count=1), 'k-')
#    axs.plot(rawspec["Wavelength"], spec1 , 'k-', lw=0.5)
#    axs.plot(rawspec["Wavelength"], spec3, 'g-', lw=0.5)
    
med = np.median(specs1, axis=0)
std = np.std(specs1, axis=0, ddof=1)
axs.set_xlabel(r'$\lambda$ ($\AA$)')
axs.set_ylabel('Arbitrary Flux')
axs.set_title("NGC 2639 Spectrum")
axs.plot(rawspec["Wavelength"], med, 'k-', lw=0.5)
axs.plot(rawspec["Wavelength"], np.median(specs1, axis=0) , 'k-', lw=0.5)
axs.fill_between(rawspec["Wavelength"], med-std, med+std, alpha=0.5, color='r', lw=0)

axs.fill_between(rawspec["Wavelength"], med-std, med+std, alpha=0.2, color='k', lw=0)
z = 0.010739
wls = [4862.721, 4960.295, 5008.239, 6549.86, 6564.614, 6585.27, 6718.29, 6732.68]
names = [r"H $\beta$", "[O III]", "[O III]", "[N II]", r"H $\alpha$", "[N II]", "[S II]", "[S II]"]
for wl, name in zip(wls, names):
    axs.axvline(wl*(1+z))
    axs.text(wl*(1+z), 0, name, fontsize=14, rotation=90,)

axs.set_xlim(6400, 7100)
axs.set_ylim(-0.0001, 0.001)




x_Halpha = []
y_Halpha = []
x_ind_Ha = np.where( (rawspec["Wavelength"] >= 6570) & (rawspec["Wavelength"] < 6700) ) 

for i in range(len(x_ind_Ha[0])):
    x_Halpha.append( rawspec["Wavelength"][x_ind_Ha[0][i]]  )
    y_Halpha.append( med[x_ind_Ha[0][i]]  ) 


Halpah, pcov = curve_fit(gaus32,
                    x_Halpha,
                    y_Halpha,
                    bounds=([0.00002,0.00002,0.0002,6605,6605,6605,1,1,1,-0.1,0.0001],[0.0004,0.0004,0.0015,6630,6630,6660,20,20,30,0.1,0.0005])  )
                    #, bounds=([], []))
                    
#a1, a2, a3, mu1, mu2, mu3, s1, s2, s3, slope, const
#
x_Ha= np.arange(x_Halpha[0],x_Halpha[-1],1) 
Ha_gua3=[]

for i in range(len(x_Ha)):
        Ha_gua3.append(gaus32(x_Ha[i],
                    Halpah[0],
                    Halpah[1],
                    Halpah[2],
                    Halpah[3],
                    Halpah[4],
                    Halpah[5],
                    Halpah[6],
                    Halpah[7],
                    Halpah[8],
                    Halpah[9],
                    Halpah[10]) )
#axs.plot(rawspec["Wavelength"], med, 'k-', lw=0.5)
#axs.plot(rawspec["Wavelength"], np.median(specs1, axis=0) , 'k-', lw=0.5)
axs.plot(x_Ha, Ha_gua3
         ,color='r', lw=2,alpha=0.7)
axs.axvline(Halpah[3],color='r', ls='--',alpha=0.5)
axs.axvline(Halpah[4],color='r', ls='--',alpha=0.5)
axs.axvline(Halpah[5],color='r', ls='--',alpha=0.5)
#
#def gauss2(x, a1, a2, mu1, mu2, s1, s2, slope, const):
#    return gauss1(x, a1, mu1, s1) + gauss1(x, a2, mu2, s2) + slope*x + const


x_SII = []
y_SII = []
x_ind_SII = np.where( (rawspec["Wavelength"] >= 6770) & (rawspec["Wavelength"] < 6820) ) 

for i in range(len(x_ind_SII[0])):
    x_SII.append( rawspec["Wavelength"][x_ind_SII[0][i]]  )
    y_SII.append( med[x_ind_SII[0][i]]  ) 


SII, pcov = curve_fit(gauss2,
                    x_SII,
                    y_SII,
                    bounds=([0.00002,0.00002,6770,6795,1,1,-0.1,0.0000],[0.00011,0.00011,6795,6820,15,10,0.1,0.0005])  )
                    #, bounds=([], [])
                    
#a1, a2, a3, mu1, mu2, mu3, s1, s2, s3, slope, const
#
x_SII= np.arange(x_SII[0],x_SII[-1],1) 
SII_gua2=[]

for i in range(len(x_SII)):
        SII_gua2.append(gauss2(x_SII[i],
                    SII[0],
                    SII[1],
                    SII[2],
                    SII[3],
                    SII[4],
                    SII[5],
                    SII[6],
                    SII[7],) )
#axs.plot(rawspec["Wavelength"], med, 'k-', lw=0.5)
#axs.plot(rawspec["Wavelength"], np.median(specs1, axis=0) , 'k-', lw=0.5)
axs.plot(x_SII, SII_gua2
         ,color='r', lw=2,alpha=0.7)
axs.axvline(SII[2],color='r', ls='--',alpha=0.5)
axs.axvline(SII[3],color='r', ls='--',alpha=0.5)


ax2.plot(rawspec["Wavelength"], senss1
         ,color='g', lw=1.5,alpha = 0.5)
#%%
Halpah




#%%
sens2 = sens_pydis(stdobsspec["Wavelength"].to_numpy(), 
               stdobsspec["Value"].to_numpy(), 
               modelspec["Wavelength"].to_numpy(), 
               modelspec["flux"].to_numpy(), 
               std_wth=np.ones_like(modelspec["Wavelength"])*5)
#%%

fig, axs = plt.subplots(1,1)



    #%%
###########################################BPT Tiral################################################
#x_log = np.arange(-2,0.3,0.03)
#x_lin = np.arange(-1,0.5,0.03) 
#y03_NII = []
#mask_y03_NII = []
#
#y01_NII = []
#
#
#yAGN_SII = []
#yLINER_SII = []
#yAGN_OI = []
#yLINER_OI = []
#
#
#y03_NII < y01_NII
#
#for i in range(len(x_log)):
#    y03_NII.append(  0.61 / ( x_log[i] -0.05 ) + 1.3 )
#    y01_NII.append(  0.61 / ( x_log[i] - 0.47) + 1.19 )
#    yAGN_SII.append( 0.72 / ( x_log[i] - 0.32) + 1.30 )   
#    yAGN_OI.append( 0.73 / ( x_log[i] + 0.59) + 1.33 )    
#    if y03_NII[i] < y01_NII[i]:
#        mask_y03_NII.append(True)
#    else: 
#        mask_y03_NII.append(False)
#        
#for i in range(len(x_lin)):
#    yLINER_OI.append( 1.89*  x_lin[i] + 0.76 )
#    yLINER_SII.append( 1.89*  x_lin[i] + 0.76 )
#    #%%
#y03_NII = np.array(y03_NII)
##%%
#x_log[mask_y03_NII]
##%%
#y03_NII[mask_y03_NII]
##%%
#np.array(y03_NII)
##%%
#x_log
##%%
#fig, axs = plt.subplots(1,1)
#axs.plot( x_log, y01_NII ,'k-', lw=1)
#axs.plot( x_log[mask_y03_NII], y03_NII[mask_y03_NII] ,'k-',ls='--', lw=1)
