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
"""Using 9 standard stars data to get sense function"""
Savepath=Path('/home/astro_02/AO2019-2/2019-10-24/Specdata')

#sds9_val = []
#sds9_wvlt = []

#SED(3814.85~5155.74 (nm)) file readiing
SEDdirpath=Path('/home/astro_02/AO2019-2/2019-10-24/SED4628')

f4628 = open(SEDdirpath/'hd4628.ctio.flx','r')
SED_4628 = f4628.readlines()

SED_wvlt = []
SED_value = []

for i in range(1193):
    wvlt = float(SED_4628[i].split(' ')[0])
    value = float(SED_4628[i].split(' ')[1])
    SED_wvlt.append(wvlt)
    SED_value.append(value)
f4628.close()


ob_val= Table.read(Savepath/'NGC676_raw.csv')['Value']
ob_wvlt = Table.read(Savepath/'NGC676_raw.csv')['Wavelength']

# (4071.71 ~ 5155.74 angstom)
#%%



#%%
sds9_val=[]
sens9 = []
spec9 = []

for i in range(9):

    sds_val = Table.read(Savepath/f'HD4628_0{i+1}.csv')['Value']
    sds_wvlt = Table.read(Savepath/f'HD4628_0{i+1}.csv')['Wavelength']
    sds9_val.append( sds_val )
    
    interp = UnivariateSpline(x=SED_wvlt, y=SED_value, s=0, k=3, ext=1)
    S_match = interp(sds_wvlt)
    sens9.append( S_match / sds_val )
    spec9.append( ob_val*sens9[i] )
#%%
fig = plt.figure(figsize=(10,5))
gs = gridspec.GridSpec(3,1)
ax = plt.subplot(gs[0:3])




ax.plot(sds_wvlt,ob_val,lw=1)

ax.grid(ls=':')
ax.legend()
plt.show()
    
    
#%%
fig = plt.figure(figsize=(10,5))
gs = gridspec.GridSpec(3,1)
ax = plt.subplot(gs[0:3])



for i in range(9):
    ax.plot(sds_wvlt,sds9_val[i],lw=1,alpha = 0.5,
            label = i
        )

ax.grid(ls=':')
ax.legend()
plt.show()
    
#%%
fig = plt.figure(figsize=(10,5))
gs = gridspec.GridSpec(3,1)
ax = plt.subplot(gs[0:3])



for i in range(9):
    ax.plot(sds_wvlt,sens9[i],lw=1,alpha = 0.5,
            label = i
        )

ax.grid(ls=':')
ax.legend()
plt.show()
#%%
fig = plt.figure(figsize=(10,5))
gs = gridspec.GridSpec(3,1)
ax = plt.subplot(gs[0:3])


for i in range(9):
    
    ax.plot(sds_wvlt,spec9[i],lw=1,alpha = 0.5,
            label = i
        )
    
ax.set_xlabel(r'$\lambda$ ($\AA$)')
ax.set_ylabel('erg/cm2/sec/Angstrom')
ax.grid(ls=':')
ax.legend()
plt.show()

#%%
# I pick the most brightest data
fig = plt.figure(figsize=(10,5))
gs = gridspec.GridSpec(3,1)
ax = plt.subplot(gs[0:3])


    
ax.plot(sds_wvlt,spec9[3],lw=1
#        ,alpha = 0.5,
#            label = i
        )
    
ax.set_xlabel(r'$\lambda$ ($\AA$)')
ax.set_ylabel('erg/cm2/sec/Angstrom')
ax.grid(ls=':')
ax.legend()
plt.show()



#%%

table=Table()

wave_length = Column(data = sds_wvlt, name='Wavelength')
Value = Column(data = spec9[3], name='Value')

table.add_column(wave_length, index = 0)
table.add_column(Value, index = 1)

table.write(Savepath/f'NGC676_spec.csv',format='ascii.csv',overwrite=True)


#%%
#####################################For HD207673 & WISEA#####################################

Savepath=Path('/home/astro_02/AO2019-2/2019-10-24/Specdata')

#sds9_val = []
#sds9_wvlt = []

#SED(3814.85~5155.74 (nm)) file readiing
SEDdirpath=Path('/home/astro_02/AO2019-2/2019-10-24/SED207673')

scc1 = CCDData.read(SEDdirpath/'hd207673_SED.fits', unit=u.adu)
scc2 = CCDData.read(SEDdirpath/'HD207673_SED2.fits', unit=u.adu)
#%%
scc1.header

#%%
SED1_wvlt = []
SED1_value = scc1.data

SED2_wvlt = []
SED2_value = scc2.data[0]


ob_val= Table.read(Savepath/'NGC2639_raw.csv')['Value']
ob_wvlt = Table.read(Savepath/'NGC2639_raw.csv')['Wavelength']
#%%#
print(len(SED2_value),len(SED1_value))
#%%
scc2.header
#%%
# I found that at 
"""
https://britastro.org/specdb/data_graph.php?obs_id=1555&obs_validated=&obs_observer_id=FJQ&r_c=1&f_c=0&o_comment=&plot=Plot
"""
for i in range(len(SED1_value)):
    SED1_wvlt.append(3934.28+ i*((7396.07-3934.28)/3905))

# And this was found at the header
for i in range(len(SED2_value)):
    SED2_wvlt.append(3500 + i*0.9)
    
ob_val= Table.read(Savepath/'WISEA_raw.csv')['Value']
ob_wvlt = Table.read(Savepath/'WISEA_raw.csv')['Wavelength']
#%%
print(SED1_wvlt[0],SED1_wvlt[-1],len(SED2_wvlt))
print(SED2_wvlt[0],SED2_wvlt[-1],len(SED2_wvlt))
#%%
sds9_val=[]
sens9 = []
spec9 = []



for i in range(9):

    
    sds_val = Table.read(Savepath/f'HD207673_0{i+1}.csv')['Value']
    sds_wvlt = Table.read(Savepath/f'HD207673_0{i+1}.csv')['Wavelength']
    sds9_val.append( sds_val )
    
    interp = UnivariateSpline(x=SED1_wvlt, y=SED1_value, s=0, k=3,ext=1)
    S_match = interp(sds_wvlt)
    sens9.append( S_match / sds_val )
    spec9.append( ob_val*sens9[i] )




#%%    
fig = plt.figure(figsize=(10,5))
gs = gridspec.GridSpec(3,1)
ax = plt.subplot(gs[0:3])
for i in range(9):
    ax.plot(sds_wvlt,sds9_val[i],lw=1,alpha = 0.5,
            label = i
        )

ax.grid(ls=':')
ax.legend()
plt.show()

#%%    
fig = plt.figure(figsize=(10,5))
gs = gridspec.GridSpec(3,1)
ax = plt.subplot(gs[0:3])
for i in range(9):
    ax.plot(sds_wvlt,sens9[i],lw=1,alpha = 0.5,
            label = i
        )

ax.grid(ls=':')
ax.legend()
plt.show()

    
#%%
fig = plt.figure(figsize=(10,5))
gs = gridspec.GridSpec(3,1)
ax = plt.subplot(gs[0:3])
for i in range(9):
    ax.plot(sds_wvlt,spec9[i],lw=1,alpha = 0.5,
            label = i
        )

ax.grid(ls=':')
ax.legend()
plt.show()
#%%


fig = plt.figure(figsize=(10,5))
gs = gridspec.GridSpec(3,1)
ax = plt.subplot(gs[0:3])

ax.plot(sds_wvlt,spec9[8],lw=1,
        )

ax.grid(ls=':')
ax.legend()
plt.show()


#%%

table=Table()

wave_length = Column(data = sds_wvlt, name='Wavelength')
Value = Column(data = spec9[8], name='Value')

table.add_column(wave_length, index = 0)
table.add_column(Value, index = 1)

table.write(Savepath/f'WISEA_spec.csv',format='ascii.csv',overwrite=True)



