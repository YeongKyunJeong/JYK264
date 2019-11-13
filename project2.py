import os
import glob

import numpy as np
from pathlib import Path
from numpy.polynomial.chebyshev import chebfit, chebval

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



#%%
newfitspath=Path('/home/astro_02/AO2019-2/2019-10-24')
ppdpath = Path(newfitspath/'PreprocessedImages')
tablepath=Path(newfitspath/'Headertable')

os.chdir(ppdpath)
allfitsname = glob.glob('*.fits')
allfitspath = list(ppdpath.glob('*.fits'))
allfitspath.sort()
allfitsname.sort()

print(len(allfitspath),"Items are searched")
#%%
tb=Table.read(tablepath/'Headersumm.csv')

#%%
print(tb.columns)
#%%
print(tb['OBJECT'])

#%%
cptb = tb[(tb['OBJECT'] == 'Comp')]
print(cptb)

#%%
for i in cptb['FILE']:
    cc = CCDData.read(i)
    cca = np.array(cc)
    print( np.where(cca == np.max(cca)) )

#%%

images = []

if os.path.exists(ppdpath/'Comp-master.fits'):
    mcomp = CCDData.read(ppdpath/'Comp-master.fits')
    print("Master comp is already exists")
    
else:
    for i in cptb['FILE']:
        cc = CCDData.read(i)
        images.append(cc)
    mcomp = combine(images, method = 'median')
    mcomp.header = cc.header
    mcomp.header.add_history(f"{len(cptb)} image(s) median combined Comparision Ramp Image)")
    mcomp = yfu.CCDData_astype(mcomp,dtype = 'float32')
    mcomp.write(ppdpath/'Comp-master.fits',overwrite=True)
#%%
mcomp.header
#%%
fig, axs = plt.subplots(1,1)
#%%
im = yfu.zimshow(axs,mcomp)


#%%
"""
fig, axs = plt.subplots(2,len(ob_list),figsize=(12,12))
for i in range(len(ob_list)):
    data1 = CCDData.read(f"{ob_list[i]}-0001.fits")
    data2 = CCDData.read(ppdpath/f"{ob_list[i]}-0001.fits")
    im1 = yfu.zimshow(axs[0][i],data1)
    im2 = yfu.zimshow(axs[1][i],data2)
    axs[0][i].set_xlabel(f"{ob_list[i]}\nraw", fontsize=8)
    axs[1][i].set_xlabel(f"{ob_list[i]}", fontsize=8)
"""
#%%
#%%
#뭔지 잘 모르겠는데 일단 추가
def disable_mplkeymaps():
    rc('keymap', 
       fullscreen='',
       home='',
       back='',
       forward='',
       pan='',
       zoom='',
       save='',
       quit='',
       grid='',
       yscale='',
       xscale='',
       all_axes=''
       )
#%%
DISPAXIS = 1 # 1 = line = python_axis_1 // 2 = column = python_axis_0
FONTSIZE = 12 # Change it on your computer if you wish.
rcParams.update({'font.size': FONTSIZE})
COMPIMAGE = os.path.join(ppdpath, 'Comp-master.fits') # Change directory if needed!
OBJIMAGE  = os.path.join(ppdpath, 'HD4628-0001.fits')
LINE_FITTER = LevMarLSQFitter()
#%%
# Parameters for IDENTIFY
FITTING_MODEL_ID = 'Chebyshev'
ORDER_ID = 4 
NSUM_ID = 10
FWHM_ID = 4 # rough guess of FWHM of lines in IDENTIFY (pixels)

# Parameters for REIDENTIFY
FITTING_MODEL_REID = 'Chebyshev' # 2-D fitting function
ORDER_SPATIAL_REID = 6
ORDER_WAVELEN_REID = 6
STEP_REID = 15  # Reidentification step size in pixels (spatial direction)
NSUM_REID = 10
TOL_REID = 5 # tolerence to lose a line in pixels

# Parameters for APALL (sky fitting and aperture extract after sky subtraction)
## parameters for finding aperture
NSUM_AP = 10
FWHM_AP = 10
STEP_AP = 10  # Recentering step size in pixels (dispersion direction)
## parameters for sky fitting
FITTING_MODEL_APSKY = 'Chebyshev'
ORDER_APSKY = 3
SIGMA_APSKY = 3
ITERS_APSKY = 5
## parameters for aperture tracing
FITTING_MODEL_APTRACE = 'Chebyshev'
ORDER_APTRACE = 3
SIGMA_APTRACE = 3
ITERS_APTRACE = 5 
# The fitting is done by SIGMA_APTRACE-sigma ITERS_APTRACE-iters clipped on the
# residual of data. 

#%%
lamphdu = fits.open(COMPIMAGE)
objhdu = fits.open(OBJIMAGE)
lampimage = lamphdu[0].data
objimage  = objhdu[0].data

if lampimage.shape != objimage.shape:
    raise ValueError('lamp and obj images should have same sizes!')

if DISPAXIS == 2:
    lampimage = lampimage.T
    objimage = objimage.T
elif DISPAXIS != 1:
    raise ValueError('DISPAXIS must be 1 or 2 (it is now {:d})'.format(DISPAXIS))

EXPTIME = objhdu[0].header['EXPTIME']
OBJNAME = objhdu[0].header['OBJECT']
# Now python axis 0 (Y-direction) is the spatial axis 
# and 1 (X-direciton) is the wavelength (dispersion) axis.
#N_SPATIAL, N_WAVELEN = np.shape(lampimage)
N_WAVELEN, N_SPATIAL = np.shape(lampimage)

N_REID = N_SPATIAL//STEP_REID # No. of reidentification
N_AP = N_WAVELEN//STEP_AP # No. of aperture finding

# ``peak_local_max`` calculates the peak location using maximum filter:
#   med1d_max = scipy.ndimage.maximum_filter(med1d, size=10, mode='constant')
# I will use this to show peaks in a primitive manner.
MINSEP_PK = 5   # minimum separation of peaks
MINAMP_PK = 0.01 # fraction of minimum amplitude (wrt maximum) to regard as peak
NMAX_PK = 50
print("setting done!")

#%%
print(np.shape(lampimage))
len(lampimage[2000])
#%%
NSUM_ID
#%%
# =============================================================================
# Identify (1): plot for manual input
# =============================================================================
# mimics IRAF IDENTIFY
#   IDENTIIFY image.fits section='middle line' nsum=NSUM_ID
lowercut_ID = N_SPATIAL//2 - NSUM_ID//2 
uppercut_ID = N_SPATIAL//2 + NSUM_ID//2
identify_1 = np.median(lampimage[:,lowercut_ID:uppercut_ID], axis=1)
"""I changed from
np.median(lampimage[lowercut_ID:uppercut_ID, :], axis=0)
        to
np.median(lampimage[:,lowercut_ID:uppercut_ID], axis=1)
"""

#%%
"""
print(lowercut_ID)
print(uppercut_ID)

len(np.median(lampimage[:,lowercut_ID:uppercut_ID], axis=1))

len(identify_1)
"""
#%%
# For plot and visualization
max_intens = np.max(identify_1)
#%%

#%%
peak_pix = peak_local_max(identify_1, indices=True, num_peaks=NMAX_PK,
                          min_distance=MINSEP_PK,
                          threshold_abs=max_intens * MINAMP_PK)
# ``peak_pix`` corresponds to the x value, since x = pixels starting from 0.

disable_mplkeymaps()
fig = plt.figure()
ax = fig.add_subplot(111)
title_str = r'Peak ($A \geq {:.2f} A_\mathrm{{max}}$, $\Delta x \geq {:.0f}$)'
# Plot original spectrum + found peak locations
x_identify = np.arange(0, len(identify_1))
ax.plot(x_identify, identify_1, lw=1)

for i in peak_pix:
    ax.plot((i, i), 
            (identify_1[i]+0.01*max_intens, identify_1[i]+0.05*max_intens),
            color='r', ls='-', lw=1)
    ax.annotate(i[0], (i, 0.9),
                xycoords = ('data', 'axes fraction'),
                fontsize='xx-small', rotation=70)
ax.grid(ls=':')
ax.set_xlabel('Pixel number')
ax.set_ylabel('Pixel value sum')
ax.set_xlim(len(identify_1), 0) # to invert x-axis
ax.set_title(title_str.format(MINAMP_PK, MINSEP_PK))
plt.show()
