import os
import glob

import numpy as np
from pathlib import Path
from numpy.polynomial.chebyshev import chebfit, chebval

from photutils import RectangularAperture
from photutils import aperture_photometry

from matplotlib import pyplot as plt
from matplotlib import gridspec, rcParams, rc
from matplotlib.widgets import Cursor

import ysfitsutilpy as yfu

import astropy.units as u
from astropy.nddata import CCDData
from astropy.table import Table, Column
from astropy.io import fits
from astropy.stats import sigma_clip, gaussian_fwhm_to_sigma
from astropy.modeling.models import Gaussian1D, Chebyshev2D
from astropy.modeling.fitting import LevMarLSQFitter

from ccdproc import Combiner, combine
from skimage.feature import peak_local_max

from scipy.interpolate import UnivariateSpline

#%%
#Now we need object/stdstars
#fitting = Xsq minimization
#interpolation(spline) pass on the all point (s = 0 (not default, should be set) -> yi(x)-S(xi) must be 0 )
#So, if there is absorption, or emission line, this affects a lot.
#Then, 1: linear spline, 2: mask that area
# Try to find SDSS spectroscopic standard star spectrum
#%%
#####################################For HD 4628 & NGC 676#####################################
# I used 9 standard stars data to get sense function and know the sense function is reliable.

Savepath=Path('/home/astro_02/AO2019-2/2019-10-24/Specdata')


# Opening SED (3814.85 ~ 5155.74 Angstrom) file readiing from flx file

SEDdirpath=Path('/home/astro_02/AO2019-2/2019-10-24/SED4628')

#f4628 = open(SEDdirpath/'hd4628.ctio.flx','r')
f4628 = open(SEDdirpath/'HD 4628 flux data.txt','r')
SED_4628 = f4628.readlines()


SED_wvlt = []
SED_value = []

for i in range(101):
    wvlt = float(SED_4628[i].split('          ')[0])
    value = float(SED_4628[i].split('          ')[1])
    adu_val = float( (SED_4628[i].split('                  ')[1]).split('\n')[0] )
    SED_wvlt.append(wvlt)
    SED_value.append( (10**(-8))*10**value)

f4628.close()

# Loading taget data(NGC 676)

ob_val= Table.read(Savepath/'NGC676_raw_01.csv')['Value']
ob_wvlt = Table.read(Savepath/'NGC676_raw_01.csv')['Wavelength']
#%%
fig = plt.figure(figsize=(10,5))
gs = gridspec.GridSpec(3,1)
ax = plt.subplot(gs[0:3])

#ax.plot(SED_wvlt,SED_adu,lw=1,alpha=0.5)
ax.plot(SED_wvlt,SED_value,lw=1)


ax.grid(ls=':')
ax.legend()
ax.set_title("Observed NGC 676")
plt.show()
#%%

sds9_val=[]
sens9 = []
spec9 = []


for i in range(9):
    # Loading standard stars data(HD 4628)
    sds_val = Table.read(Savepath/f'HD4628_0{i+1}.csv')['Value']
    sds_wvlt = Table.read(Savepath/f'HD4628_0{i+1}.csv')['Wavelength']
    sds9_val.append( sds_val )
    
    # Interpolation SED data with cubic spline, and I gave 0 to extrapolation, because that area is not reliable
    interp = UnivariateSpline(x=SED_wvlt, y=SED_value, s=0, k=3, ext=1)
    S_match = interp(sds_wvlt)
    # Sence function = Interpolated SED / obsevation data of standard star  
    # Real target intensity = sence function * observed data of target
    sens9.append( S_match / sds_val )

#%%test
    
fig = plt.figure(figsize=(10,5))
gs = gridspec.GridSpec(3,1)
ax = plt.subplot(gs[0:3])

ax.plot(sds_wvlt,S_match,lw=1)
ax.plot(SED_wvlt,SED_value,lw=1)
#ax.plot(sds_wvlt,resds,lw=1)

ax.grid(ls=':')
ax.legend()
plt.show()

#%%
# ploting observed target data, to find some emission and absorption line.
# I think we can use some way like local max finding, but not yet done.
    
fig = plt.figure(figsize=(10,5))
gs = gridspec.GridSpec(3,1)
ax = plt.subplot(gs[0:3])

ax.plot(ob_wvlt,ob_val,lw=1)

ax.grid(ls=':')
ax.legend()
ax.set_title("Observed NGC 676")
ax.set_ylabel('ADU')
plt.show()
    
#%%
# ploting all observed data of standard stars

fig = plt.figure(figsize=(10,5))
gs = gridspec.GridSpec(3,1)
ax = plt.subplot(gs[0:3])

for i in range(9):
    ax.plot(sds_wvlt,sds9_val[i],lw=1,alpha = 0.5,
            label = i
        )

ax.grid(ls=':')
ax.legend()
ax.set_title("Observed HD 4628")
ax.set_ylabel('ADU')
plt.show()
    
#%%
# ploting all sense function

fig = plt.figure(figsize=(10,5))
gs = gridspec.GridSpec(3,1)
ax = plt.subplot(gs[0:3])

for i in range(9):
    ax.plot(sds_wvlt,sens9[i],lw=1,alpha = 0.5,
            label = i
        )

ax.set_title("Sense function")
ax.set_ylabel('erg/cm2/sec/Angstrom')
ax.grid(ls=':')
ax.legend()
plt.show()
#%%















xxx=np.arange(len(sens9[0]))




#546 to 244 
ftsens=[]
ftzero_sens=[]

for i in range(9):
    ftzero_sens.append(sens9[i][0])
    ftsens.append(sens9[i][546])
#%%
print(np.where(ftsens==np.median(ftsens)), np.where(ftzero_sens==np.median(ftzero_sens)) )
#%%
fact = []
for i in range(9):
    fact.append( (sens9[1][546]-sens9[1][0])/(sens9[i][546]-sens9[i][0]) )
#%%

sense9=[]
fact
for i in range(9):
    sense9.append( (sens9[i]-sens9[i][0])*fact[i]+sens9[1][0] )
    
    for j in range(547,1101,1):
        sense9[i][j] = 0
        
#%%
fig = plt.figure(figsize=(10,5))
gs = gridspec.GridSpec(3,1)
ax = plt.subplot(gs[0:3])

for i in range(9):
    ax.plot(xxx,sense9[i],lw=1,alpha = 0.5,
            label = i
        )

ax.set_title("Sense function")
ax.set_ylabel('erg/cm2/sec/Angstrom')
ax.grid(ls=':')
ax.legend()
plt.show()
#%%

master_sens = []
for i in range(1101):
    x = []
    
    for j in range(9):
        x.append(sense9[j][i]) 
    
    master_sens.append(np.median(x))





#%%
fig = plt.figure(figsize=(10,5))
gs = gridspec.GridSpec(3,1)
ax = plt.subplot(gs[0:3])


ax.plot(xxx,master_sens,lw=1)


ax.set_title("Sense function")
ax.set_ylabel('erg/cm2/sec/Angstrom')
ax.grid(ls=':')
ax.legend()
plt.show()

#%%
table=Table()

wave_length = Column(data = sds_wvlt, name='Wavelength')
Value = Column(data = master_sens, name='Value')

table.add_column(wave_length, index = 0)
table.add_column(Value, index = 1)

table.write(Savepath/f'NGC676_MSF.csv',format='ascii.csv',overwrite=True)




#for i in range(9):
#    spec9.append( ob_val*sens9[i] )
#%%
spec = ob_val*master_sens 

#%%
# ploting all Spectrum of NGC 676

fig = plt.figure(figsize=(10,5))
gs = gridspec.GridSpec(3,1)
ax = plt.subplot(gs[0:3])

#for i in range(9):
    
#    ax.plot(sds_wvlt,spec9[i],lw=1,alpha = 0.5,
#            label = i
#        )
ax.plot(sds_wvlt,spec,lw=1) 

   
ax.set_xlabel(r'$\lambda$ ($\AA$)')
ax.set_title("Spectrum of NGC 676")
ax.set_ylabel('erg/cm2/sec/Angstrom')
ax.grid(ls=':')
ax.legend()
plt.show()


table=Table()

wave_length = Column(data = sds_wvlt, name='Wavelength')
Value = Column(data = spec, name='Value')

table.add_column(wave_length, index = 0)
table.add_column(Value, index = 1)

table.write(Savepath/f'NGC676_spectrum.csv',format='ascii.csv',overwrite=True)

#%%

fig = plt.figure(figsize=(10,5))
gs = gridspec.GridSpec(3,1)
ax = plt.subplot(gs[0:3])



spec9=[]
for i in range(9):
    ob_val= Table.read(Savepath/f'NGC676_raw_0{i+1}.csv')['Value']
    ob_wvlt = Table.read(Savepath/f'NGC676_raw_0{i+1}.csv')['Wavelength'] 
    spec9.append(ob_val*master_sens)
    
    ax.plot(ob_wvlt,spec9[i],lw=1,label = i, alpha=0.5)
    

#    table=Table()
#
#    wave_length = Column(data = ob_wvlt, name='Wavelength')
#    Value = Column(data = spec9[i], name='Value')
#
#    table.add_column(wave_length, index = 0)
#    table.add_column(Value, index = 1)
#
#    table.write(Savepath/f'NGC676_spectrum_{i+1}.csv',format='ascii.csv',overwrite=True)



ax.set_xlabel(r'$\lambda$ ($\AA$)')
ax.set_title("Spectrum of NGC 676")
ax.set_ylabel('erg/cm2/sec/Angstrom')
ax.grid(ls=':')
ax.legend()
plt.show()   
#%%
# 1st order fitting rtial
xxx=np.arange(len(spec9[0]))




#546 to 244 
ftspec=[]
ftzero_spec=[]

for i in range(9):
    ftzero_spec.append(spec9[i][0])
    ftspec.append(spec9[i][546])
#%%
print(np.where(ftspec==np.median(ftspec)),
      np.where(ftzero_spec==np.median(ftzero_spec)),
      )
#%%
fact = []
for i in range(9):
    fact.append( (spec9[7][546]-spec9[7][0])/(spec9[i][546]-spec9[i][0]) )
#%%

spect9=[]
for i in range(9):
    spect9.append( (spec9[i]-spec9[i][0])*fact[i]+spec9[1][0] )
    
    for j in range(547,1101,1):
        spect9[i][j] = 0
        
#%%
fig = plt.figure(figsize=(10,5))
gs = gridspec.GridSpec(3,1)
ax = plt.subplot(gs[0:3])

for i in range(9):
    ax.plot(xxx,spect9[i],lw=1,alpha = 0.5,
            label = i
        )

ax.set_title("Fitted Spectrum of NGC 676")
ax.set_ylabel('erg/cm2/sec/Angstrom')
ax.grid(ls=':')
ax.legend()
plt.show()
#%%

master_spec = []
for i in range(1101):
    x = []
    
    for j in range(9):
        x.append(spect9[j][i]) 
    
    master_spec.append(np.median(x))    


#%%
table=Table()

wave_length = Column(data = ob_wvlt, name='Wavelength')
Value = Column(data = master_spec, name='Value')

table.add_column(wave_length, index = 0)
table.add_column(Value, index = 1)

table.write(Savepath/f'NGC676_Mspectrum.csv',format='ascii.csv',overwrite=True)

#%%
# And next, we will try to identify emmsion and obsorption line, the intensity, and redshift of them. 

fig = plt.figure(figsize=(10,5))
gs = gridspec.GridSpec(3,1)
ax = plt.subplot(gs[0:3])

ax.plot(ob_wvlt,master_spec,lw=1)
    
ax.set_xlabel(r'$\lambda$ ($\AA$)')
ax.set_title("Spectrum of NGC676")
ax.set_ylabel('erg/cm2/sec/Angstrom')
ax.grid(ls=':')
ax.legend()
plt.show() 



#%%
# I picked the most brightest one.



#%%



#%%
#####################################For HD 207673 & WISEA#####################################
# I found the wavelength information at 
# https://britastro.org/specdb/data_graph.php?obs_id=1555&obs_validated=&obs_observer_id=FJQ&r_c=1&f_c=0&o_comment=&plot=Plot
# We couldn't find SED having unit but only relative scale, so we only can find relative flux.
# Opening SED (3934.28 ~ 7396.07 Angstrom) file readiing from fits file

Savepath=Path('/home/astro_02/AO2019-2/2019-10-24/Specdata')

SEDdirpath=Path('/home/astro_02/AO2019-2/2019-10-24/SED207673')

scc1 = CCDData.read(SEDdirpath/'hd207673_SED.fits', unit=u.adu)
#%%

SED1_wvlt = []
SED1_value = scc1.data
#%%


#%%#

for i in range(len(SED1_value)):
    SED1_wvlt.append(3934.28+ i*((7396.07-3934.28)/3905))

#%%

sds9_val=[]
sens9 = []

for i in range(9):
    sds_val = Table.read(Savepath/f'HD207673_0{i+1}.csv')['Value']
    sds_wvlt = Table.read(Savepath/f'HD207673_0{i+1}.csv')['Wavelength']
    sds9_val.append( sds_val )
    
    interp = UnivariateSpline(x=SED1_wvlt, y=SED1_value, s=0, k=3,ext=1)
    S_match = interp(sds_wvlt)
    sens9.append( (S_match / sds_val)/(np.mean(S_match[:751] / sds_val[:751])) )


#%%
fig = plt.figure(figsize=(10,5))
gs = gridspec.GridSpec(3,1)
ax = plt.subplot(gs[0:3])
for i in range(9):
    ax.plot(sds_wvlt,sds9_val[i],lw=1,alpha = 0.5,
            label = i
        )
    
ax.set_xlabel(r'$\lambda$ ($\AA$)')
ax.set_title("Observed HD 207673")
ax.set_ylabel('ADU')
ax.grid(ls=':')
ax.legend()
plt.show()

#%%    
# This is only relative intensity of sense function, so no unit.

fig = plt.figure(figsize=(10,5))
gs = gridspec.GridSpec(3,1)
ax = plt.subplot(gs[0:3])
for i in range(9):
    ax.plot(sds_wvlt,sens9[i],lw=1,alpha = 0.5,
            label = i
        )

ax.set_xlabel(r'$\lambda$ ($\AA$)')
ax.set_title("Sense function")
ax.grid(ls=':')
ax.legend()
plt.show()

#%%












xxx=np.arange(len(sens9[0]))




#546 to 244 
ftsens=[]
ftzero_sens=[]

for i in range(9):
    ftzero_sens.append(sens9[i][0])
    ftsens.append(sens9[i][750])
#%%
print(np.where(ftsens==np.median(ftsens)),
      np.where(ftzero_sens==np.median(ftzero_sens)),
      np.where(sens9[1]==0 ))
#%%
fact = []
for i in range(9):
    fact.append( (sens9[7][750]-sens9[7][0])/(sens9[i][750]-sens9[i][0]) )
#%%

sense9=[]
fact
for i in range(9):
    sense9.append( (sens9[i]-sens9[i][0])*fact[i]+sens9[1][0] )
    
    for j in range(751,1101,1):
        sense9[i][j] = 0
        
#%%
fig = plt.figure(figsize=(10,5))
gs = gridspec.GridSpec(3,1)
ax = plt.subplot(gs[0:3])

for i in range(9):
    ax.plot(xxx,sense9[i],lw=1,alpha = 0.5,
            label = i
        )

ax.set_title("Sense function")
ax.set_ylabel('Relative Intensity')
ax.grid(ls=':')
ax.legend()
plt.show()
#%%

master_sens = []
for i in range(1101):
    x = []
    
    for j in range(9):
        x.append(sense9[j][i]) 
    
    master_sens.append(np.median(x))





#%%
fig = plt.figure(figsize=(10,5))
gs = gridspec.GridSpec(3,1)
ax = plt.subplot(gs[0:3])


ax.plot(sds_wvlt,master_sens,lw=1)


ax.set_title("Sense function")
ax.set_ylabel('Relative Intensity')
ax.grid(ls=':')
ax.legend()
plt.show()

#%%
table=Table()

wave_length = Column(data = sds_wvlt, name='Wavelength')
Value = Column(data = master_sens, name='Value')

table.add_column(wave_length, index = 0)
table.add_column(Value, index = 1)

table.write(Savepath/f'WISEA_MSF.csv',format='ascii.csv',overwrite=True)



#%%
fig = plt.figure(figsize=(10,5))
gs = gridspec.GridSpec(3,1)
ax = plt.subplot(gs[0:3])

spec9=[]
for i in range(4):
    ob_val= Table.read(Savepath/f'WISEA_raw_0{i+1}.csv')['Value']
    ob_wvlt = Table.read(Savepath/f'WISEA_raw_0{i+1}.csv')['Wavelength'] 
    spec9.append(ob_val*master_sens)
    
    ax.plot(ob_wvlt,spec9[i],lw=1,label = i, alpha=0.5)
    

#    table=Table()

#    wave_length = Column(data = ob_wvlt, name='Wavelength')
#    Value = Column(data = spec9[i], name='Value')

#    table.add_column(wave_length, index = 0)
#    table.add_column(Value, index = 1)

#    table.write(Savepath/f'WISEA_spectrum_{i+1}.csv',format='ascii.csv',overwrite=True)



ax.set_xlabel(r'$\lambda$ ($\AA$)')
ax.set_title("Spectrum of WISEA")
ax.set_ylabel('ADU')
ax.grid(ls=':')
ax.legend()
plt.show()   
    
#%%
fig = plt.figure(figsize=(10,5))
gs = gridspec.GridSpec(3,1)
ax = plt.subplot(gs[0:3])


for i in range(4):
    ob_val= Table.read(Savepath/f'WISEA_raw_0{i+1}.csv')['Value']
    ob_wvlt = Table.read(Savepath/f'WISEA_raw_0{i+1}.csv')['Wavelength'] 
    
    ax.plot(ob_wvlt,ob_val,lw=1,label = i,alpha=0.5)
    
ax.set_xlabel(r'$\lambda$ ($\AA$)')
ax.set_title("Spectrum of WISEA")
ax.set_ylabel('ADU')
ax.grid(ls=':')
ax.legend()
plt.show() 
#%%

master_spec = []
for i in range(1101):
    x = []
    
    for j in range(4):
        x.append(spec9[j][i]) 
    
    master_spec.append(np.median(x))
#%%
table=Table()

wave_length = Column(data = ob_wvlt, name='Wavelength')
Value = Column(data = master_spec, name='Value')

table.add_column(wave_length, index = 0)
table.add_column(Value, index = 1)

table.write(Savepath/f'WISEA_Mspectrum.csv',format='ascii.csv',overwrite=True)

#%%
# And next, we will try to identify emmsion and obsorption line, the intensity, and redshift of them. 

fig = plt.figure(figsize=(10,5))
gs = gridspec.GridSpec(3,1)
ax = plt.subplot(gs[0:3])

ax.plot(ob_wvlt,master_spec,lw=1)
    
ax.set_xlabel(r'$\lambda$ ($\AA$)')
ax.set_title("Spectrum of WISEA")
ax.set_ylabel('Relative Intensity')
ax.grid(ls=':')
ax.legend()
plt.show() 



