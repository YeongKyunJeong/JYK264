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
DISPAXIS = 2 # 1 = line = python_axis_1 // 2 = column = python_axis_0
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
N_SPATIAL, N_WAVELEN = np.shape(lampimage)

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
# =============================================================================
# Identify (1): plot for manual input
# =============================================================================
# mimics IRAF IDENTIFY
#   IDENTIIFY image.fits section='middle line' nsum=NSUM_ID
lowercut_ID = N_SPATIAL//2 - NSUM_ID//2 
uppercut_ID = N_SPATIAL//2 + NSUM_ID//2
identify_1 = np.median(lampimage[lowercut_ID:uppercut_ID, :], axis=0)

# For plot and visualization
max_intens = np.max(identify_1)

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

#%%

ID_init = dict(pixel_init=[1201, 1180,
                           1160, 1134, 
                           1117, 1106,
                           1091,
                           1051, 1030, 1012,
                           1003, 954, 931, 898, 882, 838],
    
               wavelength=[5852.49, 5944.83,
                           6030.00, 6143.06,
                           6217.28, 6266.49,
                           6334.43,
                           6506.53, 6598.95, 6678.28,
                           6717.01, 6929.47, 7032.41, 7173.94, 7245.17, 7438.90])

#%%
print(len(ID_init['pixel_init'])
        ,len(ID_init['wavelength']))
#%%
# Make it as astropy table
ID_init = Table(ID_init)

peak_gauss = []
fitter = LevMarLSQFitter()

for peak_pix in ID_init['pixel_init']:
    #TODO: put something like "lost_factor" to multiply to FWHM_ID in the bounds.
    g_init = Gaussian1D(amplitude = identify_1[peak_pix], 
                       mean = peak_pix, 
                       stddev = FWHM_ID * gaussian_fwhm_to_sigma,
                       bounds={'amplitude': (0, 2*identify_1[peak_pix]),
                               'mean':(peak_pix-FWHM_ID, peak_pix+FWHM_ID),
                               'stddev':(0, FWHM_ID)})
    fitted = LINE_FITTER(g_init, x_identify, identify_1)
    peak_gauss.append(fitted.mean.value)
    
peak_gauss = Column(data=peak_gauss, name='pixel_gauss', dtype=float)
peak_shift = Column(data=peak_gauss - ID_init['pixel_init'],
                    name='pixel_shift', dtype=float)
ID_init["pixel_gauss"] = peak_gauss
ID_init["pixel_shift"] = peak_gauss - ID_init['pixel_init']
ID_init.sort('wavelength')
ID_init.pprint()
# If you want to save:
#ID_init.write('user_input_identify.csv', 
#              format='ascii.csv', overwrite=True)
#%%
if FITTING_MODEL_ID.lower() == 'chebyshev':
    coeff_ID, fitfull = chebfit(ID_init['pixel_gauss'], 
                                ID_init['wavelength'], 
                                deg=ORDER_ID,
                                full=True)
    fitRMS = np.sqrt(fitfull[0][0]/len(ID_init))
    rough_error = ( np.ptp(ID_init['wavelength']) 
                   / np.ptp(ID_init['pixel_gauss']) ) / 2
    # rough_error = (wave_max - wave_min) / (spatial_max - spatial_min)
    residual = ( ID_init['wavelength'] 
                - chebval(ID_init['pixel_gauss'], coeff_ID))
    res_range = np.max(np.abs(residual))
else:
    raise ValueError('Function {:s} is not implemented.'.format(FITTING_MODEL_REID))

fig = plt.figure(figsize=(10, 5))
gs = gridspec.GridSpec(3, 1)
ax1 = plt.subplot(gs[0:2])
ax2 = plt.subplot(gs[2])
plt.setp(ax1.get_xticklabels(), visible=False)
ax1.plot(ID_init['pixel_gauss'], ID_init['wavelength'],
         ls=':', color='k', ms=10, marker='+')
ax2.plot(ID_init['pixel_gauss'], residual, 
         ls='', color='k', ms=10, marker='+')
ax2.axhline(y=0, color='b', ls=':')
# rough error ~ +- wavelength resolution/2
ax2.axhline(y=-rough_error/2, color='r', ls=':')
ax2.axhline(y=+rough_error/2, color='r', ls=':')
ax2.set_ylim(min(-rough_error/2 * 1.1, -1.1*res_range), 
             max(+rough_error/2 * 1.1, +1.1*res_range))
ax1.set_ylabel(r'Wavelength ($\AA$)')
ax2.set_ylabel(r'Residual ($\AA$)')
ax2.set_xlabel('Pixel along dispersion axis')
ax1.grid(ls=':')
ax2.grid(ls=':')
plt.suptitle(('First Identify (Chebyshev order {:d})\n'.format(ORDER_ID) 
              + r'RMSE = {:.2f} $\AA$'.format(fitRMS)))
plt.show()
#%%
# =============================================================================
# reidentify
# =============================================================================
# For REIDENTIFY, I used astropy.modeling.models.Chebyshev2D

line_REID = np.zeros((len(ID_init), N_REID-1))
spatialcoord = np.arange(0, (N_REID - 1) * STEP_REID, STEP_REID) + STEP_REID / 2

print('Reidentify each section by Chebyshev (order {:d})'.format(ORDER_ID))
print('section      |  found  |  RMS')

for i in range(0, N_REID-1):
    lower_cut, upper_cut = i*STEP_REID, (i+1)*STEP_REID
    reidentify_i = np.sum(lampimage[lower_cut:upper_cut, :], 
                          axis=0)
    peak_gauss_REID = []
    
    for peak_pix_init in ID_init['pixel_gauss']:
        # TODO: Will a simple astropy fitting work well for general cases?
        #TODO: put something like "lost_factor" to multiply to FWHM_ID in the bounds.
        search_min = int(np.around(peak_pix_init - TOL_REID))
        search_max = int(np.around(peak_pix_init + TOL_REID))
        cropped = reidentify_i[search_min:search_max]
        x_cropped = np.arange(len(cropped)) + search_min
        
        #TODO: put something like "lost_factor" to multiply to FWHM_ID in the bounds.
        A_init = np.max(cropped)
        mean_init = peak_pix_init
        stddev_init = FWHM_ID * gaussian_fwhm_to_sigma
        g_init = Gaussian1D(amplitude=A_init, mean=mean_init, stddev=stddev_init,
                            bounds={'amplitude':(0, 2*np.max(cropped)) ,
                                    'stddev':(0, TOL_REID)})
        g_fit = fitter(g_init, x_cropped, cropped)    
        fit_center = g_fit.mean.value
        if abs(fit_center - peak_pix_init) > TOL_REID:
            peak_gauss_REID.append(np.nan)
            continue
        peak_gauss_REID.append(fit_center)

    peak_gauss_REID = np.array(peak_gauss_REID)
    nonan_REID = np.isfinite(peak_gauss_REID)
    line_REID[:, i] = peak_gauss_REID
    peak_gauss_REID_nonan = peak_gauss_REID[nonan_REID]
    n_tot = len(peak_gauss_REID)
    n_found = np.count_nonzero(nonan_REID)
    
    if FITTING_MODEL_ID.lower() == 'chebyshev':
        coeff_REID1D, fitfull = chebfit(peak_gauss_REID_nonan,
                                        ID_init['wavelength'][nonan_REID], 
                                        deg=ORDER_WAVELEN_REID,
                                        full=True)
        fitRMS = np.sqrt(fitfull[0][0]/n_found)
    
    else:
        raise ValueError('Function {:s} is not implemented.'.format(FITTING_MODEL_REID))

    print('[{:04d}:{:04d}]\t{:d}/{:d}\t{:.3f}'.format(lower_cut, upper_cut,
                                                      n_found, n_tot, fitRMS))
#%%
    