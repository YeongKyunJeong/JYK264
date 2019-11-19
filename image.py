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

images = []


#%%
DISPAXIS = 2 # 1 = line = python_axis_1 // 2 = column = python_axis_0
FONTSIZE = 12 # Change it on your computer if you wish.
rcParams.update({'font.size': FONTSIZE})
COMPIMAGE = os.path.join(ppdpath, 'Comp-master.fits') # Change directory if needed!
OBJIMAGE  = os.path.join(newfitspath, 'HD207673-0001.fits')
LINE_FITTER = LevMarLSQFitter()

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
GAIN = objhdu[0].header['EGAIN']

# ``peak_local_max`` calculates the peak location using maximum filter:
#   med1d_max = scipy.ndimage.maximum_filter(med1d, size=10, mode='constant')
# I will use this to show peaks in a primitive manner.
MINSEP_PK = 5   # minimum separation of peaks
MINAMP_PK = 0.01 # fraction of minimum amplitude (wrt maximum) to regard as peak
NMAX_PK = 50
print("setting done!")
#%%
fig, axs = plt.subplots(1,1)
im = yfu.zimshow(axs,objimage)
