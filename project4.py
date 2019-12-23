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
Savepath=Path('/home/astro_02/AO2019-2/2019-10-24/Specdata')
# Opening SED (3814.85 ~ 5155.74 Angstrom) file readiing from flx file

SEDdirpath=Path('/home/astro_02/AO2019-2/2019-10-24/SED207673')


# ``peak_local_max`` calculates the peak location using maximum filter:
#   med1d_max = scipy.ndimage.maximum_filter(med1d, size=10, mode='constant')
# I will use this to show peaks in a primitive manner.
MINSEP_PK = 4   # minimum separation of peaks
MINAMP_PK = 0.001 # fraction of minimum amplitude (wrt maximum) to regard as peak
NMAX_PK = 30

FITTING_MODEL_ID = 'Chebyshev'
ORDER_ID = 4 
NSUM_ID = 10
FWHM_ID = 4 

raw_val= Table.read(Savepath/'WISEA_raw.csv')['Value']
raw_wvlt = Table.read(Savepath/'WISEA_raw.csv')['Wavelength']

ob_val= Table.read(Savepath/'WISEA_spec(rel).csv')['Value']
ob_wvlt = Table.read(Savepath/'WISEA_spec(rel).csv')['Wavelength']

print("setting done!")
#%%
# =============================================================================
# Identify (1): plot for manual input
# =============================================================================
# mimics IRAF IDENTIFY
#   IDENTIIFY image.fits section='middle line' nsum=NSUM_ID
#%%
zeromask = ( (ob_val) ==0 )
identify_1 = (ob_val[~zeromask])
max_intens = np.max(identify_1)
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
x_identify = np.arange(len(identify_1))
ax.plot(x_identify, identify_1, lw=1)

for i in peak_pix:
    ax.plot((i, i), 
            (identify_1[i]+0.01*max_intens, identify_1[i]+0.05*max_intens),
            color='r', ls='-', lw=1)
    ax.annotate(i[0], (i, 0.9),
                xycoords = ('data', 'axes fraction'),
                fontsize='xx-small', rotation=70)
ax.grid(ls=':')
ax.set_xlabel('wavelenth')
ax.set_ylabel('rel intensity')
#ax.set_xlim(len(identify_1), 0) # to invert x-axis
ax.set_title(title_str.format(MINAMP_PK, MINSEP_PK))
plt.show()



#%%
"""Image presentation"""
fig = plt.figure(figsize=(10,5))
gs = gridspec.GridSpec(3,1)
ax = plt.subplot(gs[0:3])

ax.plot(raw_wvlt,raw_val,lw=1,
        )

ax.set_xlabel(r'$\lambda$ ($\AA$)')
ax.set_title("Relative Intensity Spectrum of WISEA")
ax.grid(ls=':')
ax.legend()
plt.show()
