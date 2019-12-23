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
newfitspath=Path('/home/astro_02/AO2019-2/2019-10-24')
ppdpath = Path(newfitspath/'PreprocessedImages')
tablepath=Path(newfitspath/'Headertable')

os.chdir(ppdpath)
allfitsname = glob.glob('*.fits')
allfitspath = list(ppdpath.glob('*.fits'))
allfitspath.sort()
allfitsname.sort()

SKYSBPATH=Path('/home/astro_02/AO2019-2/2019-10-24/SkysubImages')

print(len(allfitspath),"Items are searched")
#%%
tb=Table.read(tablepath/'Headersumm.csv')

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
COMPIMAGE = ppdpath/'Comp-master.fits' # Change directory if needed!
OBJIMAGE  = ppdpath/'NGC676-0001.fits'
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
#for NGC2639
#ORDER_APTRACE = 9
#SIGMA_APTRACE = 9
#ITERS_APTRACE = 10

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
#fig, axs = plt.subplots(1,1)
#im = yfu.zimshow(axs,objimage)
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


#ax.plot(x_identify[10:850], identify_1[10:850], lw=1)
#
#pixel_list = [804, 702, 681,
#                           661, 635, 
#                           618, 607,
#                           592,
#                           552, 531, 513,
#                           504, 455, 432, 399, 383, 339,
#                           127, 100, 36]

    
#for i in pixel_list:
#    ax.plot((i, i), 
#            (identify_1[i]+0.01*max_intens, identify_1[i]+0.05*max_intens),
#            color='r', ls='-', lw=1)
#    ax.annotate(i, (i, 0.9),
#                xycoords = ('data', 'axes fraction'),
#                fontsize = 15,
#                rotation=70)
    
ax.grid(ls=':')
ax.set_xlabel('Pixel number')
ax.set_ylabel('Pixel value sum')
#ax.set_xlim(len(identify_1[10:850]), 0) # to invert x-axis
ax.set_xlim(len(identify_1), 0)
ax.set_title(title_str.format(MINAMP_PK, MINSEP_PK))
plt.show()
#%%
ID_init = dict(pixel_init=[804, 702, 681,
                           661, 635, 
                           618, 607,
                           592,
                           552, 531, 513,
                           504, 455, 432, 399, 383, 339,
                           127, 100, 36],
    
               wavelength=[5400.5620,5852.49, 5944.83,
                           6030.00, 6143.06,
                           6217.28, 6266.49,
                           6334.43,
                           6506.53, 6598.95, 6678.28,
                           6717.01, 6929.47, 7032.41, 7173.94, 7245.17, 7438.90,
                           8377.6070,8495.3600,8780.6220])

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
points = np.vstack((line_REID.flatten(),
                    np.tile(spatialcoord, len(line_REID))))

points = points.T # list of ()

nanmask = ( np.isnan(points[:,0]) | np.isnan(points[:,1]) )
points = points[~nanmask]
values = np.repeat(ID_init['wavelength'], N_REID - 1)
values = np.array(values.tolist())
values = values[~nanmask]
errors = np.ones_like(values)

if FITTING_MODEL_REID.lower() == 'chebyshev':
    coeff_init = Chebyshev2D(x_degree=ORDER_WAVELEN_REID, y_degree=ORDER_SPATIAL_REID)
    fit2D_REID = fitter(coeff_init, points[:, 0], points[:, 1], values)
    ww, ss = np.mgrid[:N_WAVELEN, :N_SPATIAL]
else:
    raise ValueError('Function {:s} is not implemented.'.format(FITTING_MODEL_REID))
    

#%%
fig = plt.figure(figsize=(10, 6))
gs = gridspec.GridSpec(3, 3)
ax1 = plt.subplot(gs[:2, :2])
ax2 = plt.subplot(gs[2, :2])
ax3 = plt.subplot(gs[:2, 2])
#plt.setp(ax2.get_xticklabels(), visible=False)
#plt.setp(ax3.get_yticklabels(), visible=False)

title_str = ('Reidentify and Wavelength Map\n'
             + 'func=Chebyshev, order (wavelength, dispersion) = ({:d}, {:d})')
plt.suptitle(title_str.format(ORDER_WAVELEN_REID, ORDER_SPATIAL_REID))

interp_min = line_REID[~np.isnan(line_REID)].min()
interp_max = line_REID[~np.isnan(line_REID)].max()

ax1.imshow(fit2D_REID(ww, ss).T, origin='lower')
ax1.axvline(interp_max, color='r', lw=1)
ax1.axvline(interp_min, color='r', lw=1)
#ax1.text((interp_max+interp_min)/2, -10, "Fitting Region",
#         horizontalalignment="center",
#         bbox={'facecolor':'red', 'alpha':0.8, 'pad':10})
ax1.plot(points[:, 0], points[:, 1], ls='', marker='+', color='r',
         alpha=0.8, ms=10)


for i in (1, 2, 3):
    vcut = N_WAVELEN * i/4
    hcut = N_SPATIAL * i/4
    vcutax  = np.arange(0, N_SPATIAL, STEP_REID) + STEP_REID/2
    hcutax  = np.arange(0, N_WAVELEN, 1)
    vcutrep = np.repeat(vcut, len(vcutax))
    hcutrep = np.repeat(hcut, len(hcutax))
    
    ax1.axvline(x=vcut, ls=':', color='k')   
    ax1.axhline(y=hcut, ls=':', color='k')
    
    ax2.plot(hcutax, fit2D_REID(hcutax, hcutrep), lw=1, 
             label="hcut {:d}".format(int(hcut)))

    vcut_profile = fit2D_REID(vcutrep, vcutax)
    vcut_normalize = vcut_profile / np.median(vcut_profile)
    ax3.plot(vcut_normalize, vcutax, lw=1,
             label="vcut {:d}".format(int(vcut)))

ax1.set_ylabel('Spatial direction')
ax2.grid(ls=':')
ax2.legend()
ax2.set_xlabel('Dispersion direction')
ax2.set_ylabel('Wavelength change\n(horizontal cut)')

ax3.axvline(1, ls=':', color='k')
ax3.grid(ls=':', which='both')
ax3.set_xlabel('Fractional change (vertical cut)')
ax3.legend()

ax1.set_ylim(0, N_SPATIAL)
ax1.set_xlim(0, N_WAVELEN)
ax2.set_xlim(0, N_WAVELEN)
ax3.set_ylim(0, N_SPATIAL)
plt.show()

#%%
# =============================================================================
# apall (1): Plot a cut
# =============================================================================
lower_cut = N_WAVELEN//2 - NSUM_AP//2 
upper_cut = N_WAVELEN//2 + NSUM_AP//2
apall_1 = np.sum(objimage[:, lower_cut:upper_cut], axis=1)
max_intens = np.max(apall_1)

peak_pix = peak_local_max(apall_1, indices=True, num_peaks=10,
                          min_distance=int(N_SPATIAL/100),
                          threshold_abs=np.median(apall_1))

x_apall = np.arange(0, len(apall_1))


fig = plt.figure()
ax = fig.add_subplot(111)
title_str = r'Peak ($A \geq {:.2f} A_\mathrm{{max}}$, $\Delta x \geq {:.0f}$)'

ax.plot(x_apall, apall_1, lw=1)

for i in peak_pix:
    ax.plot((i, i), 
            (apall_1[i]+0.01*max_intens, apall_1[i]+0.05*max_intens),
            color='r', ls='-', lw=1)
    ax.annotate(i[0], (i, 0.9),
                xycoords = ('data', 'axes fraction'),
                fontsize='xx-small', rotation=70)
ax.grid(ls=':')
ax.set_xlabel('Pixel number')
ax.set_ylabel('Pixel value')
ax.set_xlim(0, len(apall_1))
ax.set_title(title_str.format(np.median(apall_1), int(N_SPATIAL/100)))
plt.show()

#%%
# =============================================================================
# apall(2): manually select sky, see how the fit works
# =============================================================================
ap_init = int(peak_pix[0])

#for original sky region(NGC676, HD207673, HD4628)
ap_sky = np.array([ap_init-40, ap_init-20, ap_init+20, ap_init+40])
#for NGC 2639
#ap_sky = np.array([ap_init-55, ap_init-35, ap_init+20, ap_init+40])
#for WISEA
#ap_sky = np.array([ap_init-40, ap_init-20, ap_init+35, ap_init+55])

# Regions to use as sky background. xl1 - 1, xu1, xl2 - 1, xu2. (0-indexing)
#   Sky region should also move with aperture center!
#   from ``ap_center - 50`` to ``ap_center - 40``, for example, should be used.
# TODO: accept any even-number-sized ap_sky.

# Interactive check
x_sky = np.hstack( (np.arange(ap_sky[0], ap_sky[1]), 
                    np.arange(ap_sky[2], ap_sky[3])))
sky_val = np.hstack( (apall_1[ap_sky[0]:ap_sky[1]], 
                      apall_1[ap_sky[2]:ap_sky[3]]))


fig = plt.figure(figsize=(8, 5))
ax = fig.add_subplot(111)
title_str = r'Skyfit: {:s} order {:d} ({:.1f}-sigma {:d}-iters)'
ax.plot(x_apall, apall_1, lw=1)

if FITTING_MODEL_APSKY.lower() == 'chebyshev':
    # TODO: maybe the clip is "3-sigma clip to residual and re-fit many times"?
    clip_mask = sigma_clip(sky_val, sigma=SIGMA_APSKY, iters=ITERS_APSKY).mask
    coeff_apsky, fitfull = chebfit(x_sky[~clip_mask], 
                                   sky_val[~clip_mask],
                                   deg=ORDER_APSKY,
                                   full=True)
    fitRMS = np.sqrt(fitfull[0][0]/n_found)
    n_sky = len(x_sky)
    n_rej = np.count_nonzero(clip_mask)
    sky_fit = chebval(x_apall, coeff_apsky) 
    ax.plot(x_apall, sky_fit, ls='--',
            label='Sky Fit ({:d}/{:d} used)'.format(n_sky - n_rej, n_sky))
    ax.plot(x_sky[clip_mask], sky_val[clip_mask], marker='x', ls='', ms=10)
    [ax.axvline(i, lw=1, color='k', ls='--') for i in ap_sky]
    ax.legend()
    ax.set_title(title_str.format(FITTING_MODEL_APSKY, ORDER_APSKY,
                                  SIGMA_APSKY, ITERS_APSKY))
    ax.set_xlabel('Pixel number')
    ax.set_ylabel('Pixel value sum')
    
else:
    raise ValueError('Function {:s} is not implemented.'.format(FITTING_MODEL_REID))
ax.grid(ls=':')
ax.set_xlabel('Pixel number')
ax.set_ylabel('Pixel value')
plt.show()
#%%
#xx= np.arange(221)
#
#fig = plt.figure(figsize=(8, 5))
#ax = fig.add_subplot(111)
#ax.plot(xx, sky_fit, lw=1)
#%%

# =============================================================================
# apall (3): aperture trace
# =============================================================================
# within +- 100 pixels around the aperture, the wavelength does not change much
# as can be seen from reidentify figure 
# (in NHAO case, ~ 0.01% ~ 0.1 Angstrom order).
# So it's safe to assume the wavelength is constant over around such region,
# (spatial direction) and thus I will do sky fitting from this column,
# without considering the wavelength change along a column.
# Then aperture extraction will map the pixel to wavelength using aperture
# trace solution.

aptrace = []
aptrace_fwhm = []
#coeff_apsky = []
#aptrace_apsum = []
#aptrace_wavelen = []
# TODO: This is quite slow as for loop used: improvement needed.
# I guess the problem is sigma-clipping rather than fitting process..
for i in range(N_AP - 1):
    lower_cut, upper_cut = i*STEP_AP, (i+1)*STEP_AP
    
    apall_i = np.sum(objimage[:, lower_cut:upper_cut], axis=1)
    sky_val = np.hstack( (apall_i[ap_sky[0]:ap_sky[1]], 
                          apall_i[ap_sky[2]:ap_sky[3]]))

    # Subtract fitted sky
    if FITTING_MODEL_APSKY.lower() == 'chebyshev':
        # TODO: maybe we can put smoothing function as IRAF APALL's b_naverage 
        clip_mask = sigma_clip(sky_val, sigma=SIGMA_APSKY, iters=ITERS_APSKY).mask
        coeff, fitfull = chebfit(x_sky[~clip_mask], 
                                 sky_val[~clip_mask],
                                 deg=ORDER_APSKY,
                                 full=True)
        apall_i -= chebval(x_apall, coeff)
#        fitRMS = np.sqrt(fitfull[0][0]/n_found)
#        n_sky = len(x_sky)
#        n_rej = np.count_nonzero(clip_mask)
    
    else:
        raise ValueError('Function {:s} is not implemented.'.format(FITTING_MODEL_APSKY))

    #TODO: put something like "lost_factor" to multiply to FWHM_ID in the bounds.
    search_min = int(np.around(ap_init - 3*FWHM_AP))
    search_max = int(np.around(ap_init + 3*FWHM_AP))
    cropped = apall_i[search_min:search_max]
    x_cropped = np.arange(len(cropped))
    peak_pix = peak_local_max(cropped, 
                              min_distance=FWHM_AP,
                              indices=True,
                              num_peaks=1)
    if len(peak_pix) == 0:
        aptrace.append(np.nan)
        continue
    peak_pix = peak_pix[0][0]
    
    #TODO: put something like "lost_factor" to multiply to FWHM_ID in the bounds.
    g_init = Gaussian1D(amplitude=cropped[peak_pix], 
                       mean=peak_pix, 
                       stddev=FWHM_AP * gaussian_fwhm_to_sigma,
                       bounds={'amplitude':(0, 2*cropped[peak_pix]) ,
                               'mean':(peak_pix-3*FWHM_AP, peak_pix+3*FWHM_AP),
                               'stddev':(0, FWHM_AP)})
    fitted = fitter(g_init, x_cropped, cropped)
    center_pix = fitted.mean.value + search_min
    std_pix = fitted.stddev.value
    aptrace_fwhm.append(fitted.fwhm)
    aptrace.append(center_pix)
#    coeff_apsky.append(coeff)
#    aptrace_apsum.append(apsum)
#    apsum_lower = int(np.around(center_pix - apsum_sigma_lower * std_pix))
#    apsum_upper = int(np.around(center_pix + apsum_sigma_upper * std_pix))
#    apsum = np.sum(apall_i[apsum_lower:apsum_upper])

aptrace = np.array(aptrace)
aptrace_fwhm = np.array(aptrace_fwhm)
#coeff_apsky = np.array(coeff_apsky)
#aptrace_apsum = np.array(aptrace_apsum)

#%%
# =============================================================================
# apall(4): aperture trace fit
# =============================================================================
mask = np.isnan(aptrace)

x_aptrace = np.arange(N_AP-1) * STEP_AP
coeff_aptrace = chebfit(x_aptrace[~mask], aptrace[~mask], deg=ORDER_APTRACE)

#aptrace_finite = aptrace[np.isfinite(aptrace)]
#aptrace_finite

resid_mask = sigma_clip(aptrace - chebval(x_aptrace, coeff_aptrace), 
                        sigma=SIGMA_APTRACE, iters=ITERS_APTRACE).mask

x_aptrace_fin = x_aptrace[~resid_mask]
aptrace_fin = aptrace[~resid_mask]
coeff_aptrace_fin = chebfit(x_aptrace_fin, aptrace_fin, deg=ORDER_APTRACE)
fit_aptrace_fin   = chebval(x_aptrace_fin, coeff_aptrace_fin)
resid_aptrace_fin = aptrace_fin - fit_aptrace_fin
del_aptrace = ~np.in1d(x_aptrace, x_aptrace_fin) # deleted points

fig = plt.figure(figsize=(10,8))
gs = gridspec.GridSpec(3, 1)
ax1 = plt.subplot(gs[0:2], sharex=ax2)
ax2 = plt.subplot(gs[2])

title_str = ('Aperture Trace Fit ({:s} order {:d})\n'
            + 'Residuials {:.1f}-sigma, {:d}-iters clipped')
plt.suptitle(title_str.format(FITTING_MODEL_APTRACE, ORDER_APTRACE,
                              SIGMA_APTRACE, ITERS_APTRACE))
ax1.plot(x_aptrace, aptrace, ls='', marker='+', ms=10)
ax1.plot(x_aptrace_fin, fit_aptrace_fin, ls='--',
         label="Aperture Trace ({:d}/{:d} used)".format(len(aptrace_fin), N_AP-1))
ax1.plot(x_aptrace[del_aptrace], aptrace[del_aptrace], ls='', marker='x', ms=10)
ax1.legend()
ax2.plot(x_aptrace_fin, resid_aptrace_fin, ls='', marker='+')
#ax2.plot(x_aptrace, aptrace - chebval(x_aptrace, coeff_aptrace_fin), 
#         ls='', marker='+')
ax2.axhline(+np.std(resid_aptrace_fin, ddof=1), ls=':', color='k')
ax2.axhline(-np.std(resid_aptrace_fin, ddof=1), ls=':', color='k', 
            label='residual std')

ax1.set_ylabel('Found object position')
ax2.set_ylabel('Residual (pixel)')
ax2.set_xlabel('Dispersion axis (pixel)')
ax1.grid(ls=':')
ax2.grid(ls=':')
ax2.set_ylim(-.5, .5)
ax2.legend()
plt.show()
#plt.savefig('aptrace.png', bbox_inches='tight')

#%%# =============================================================================
# apall(5): aperture sum (NEW)
# =============================================================================
#for standard star(HD4628, HD61931, HD207673, +NGC676)
apsum_sigma_lower = 5
apsum_sigma_upper = 5
#for object(extended source: NGC2639, WISEA)
#apsum_sigma_lower = 3 # See below
#apsum_sigma_upper = 3 
# lower and upper limits of aperture to set from the center in gauss-sigma unit.
ap_fwhm = np.median(aptrace_fwhm[~resid_mask])
ap_sigma = ap_fwhm * gaussian_fwhm_to_sigma
apheight = (apsum_sigma_lower + apsum_sigma_upper)*ap_sigma

x_ap = np.arange(N_WAVELEN)
y_ap = chebval(x_ap, coeff_aptrace_fin)
ap_wavelen = fit2D_REID(x_ap, y_ap)
ap_sky_offset = ap_sky - ap_init

ap_summed  = []
data_skysub = []

for i in range(N_WAVELEN):
    cut_i = objimage[:, i].copy()
    ap_sky_i = int(y_ap[i]) + ap_sky_offset
   
    x_sky = np.hstack( (np.arange(ap_sky_i[0], ap_sky_i[1]),
                        np.arange(ap_sky_i[2], ap_sky_i[3])))
    sky_val = np.hstack( (cut_i[ap_sky_i[0]:ap_sky_i[1]],
                          cut_i[ap_sky_i[2]:ap_sky_i[3]]))
    clip_mask = sigma_clip(sky_val, sigma=SIGMA_APSKY, maxiters=ITERS_APSKY).mask
   
    coeff = chebfit(x_sky[~clip_mask],
                    sky_val[~clip_mask],
                    deg=ORDER_APSKY)

    data_skysub.append(cut_i - chebval(np.arange(cut_i.shape[0]), coeff))

data_skysub = np.array(data_skysub).T
hdr = objhdu[0].header
hdr.add_history(f"Sky subtracted using sky offset = {ap_sky_offset}, "
                + f"{SIGMA_APSKY}-sigma {ITERS_APSKY}-iter clipping "
                + f"to fit order {ORDER_APSKY} Chebyshev")
_ = fits.PrimaryHDU(data=data_skysub, header=hdr)
_.data = _.data.astype('float32')
_.writeto(SKYSBPATH/(OBJIMAGE.stem+".skysub.fits"), overwrite=True)

pos = np.array([x_ap, y_ap]).T
aps = RectangularAperture(positions=pos, w=1, h=apheight, theta=0)
phot = aperture_photometry(data_skysub, aps, method='subpixel', subpixels=30)
ap_summed = phot['aperture_sum'] / EXPTIME

fig, axs = plt.subplots(2, 1, figsize=(10, 6), sharex=False, sharey=False, gridspec_kw=None)

axs[0].imshow(objimage, vmin=0, vmax=4900, origin='lower')
axs[1].imshow(data_skysub, vmin=0, vmax=200, origin='lower')
axs[1].plot(pos[:,0], pos[:,1], 'r-')
axs[0].set(title="Before sky subtraction")
axs[1].set(title="Sky subtracted")

for ax in [axs[0], axs[1]]:
    ax.plot(x_ap, y_ap+apsum_sigma_upper*ap_sigma, 'r-', lw=1)
    ax.plot(x_ap, y_ap-apsum_sigma_lower*ap_sigma, 'r-', lw=1)
#    ax.set(ylim=(400, 700))
    for offset in ap_sky_offset:
        ax.plot(x_ap, y_ap+offset, 'y-', lw=1)
    
plt.tight_layout()
plt.show()



#%%

fig = plt.figure(figsize=(10,5))
gs = gridspec.GridSpec(3,1)
ax = plt.subplot(gs[0:3])

#plt.setp(ax.get_xticklabels(), visible=False)
title_str = ('{:s}\nobject = {:s}, EXPTIME = {:.1f} s')
label_str = r'aperture width $\approx$ {:.1f} pix'
plt.suptitle(title_str.format(OBJIMAGE.name,
                              OBJNAME, EXPTIME))
ax.plot(ap_wavelen, ap_summed, lw = 1,
         label=label_str.format(apheight))
ax.set_ylabel('Instrument Intensity\n(apsum/EXPTIME)')
ax.set_xlabel(r'$\lambda$ ($\AA$)')

ax.legend()
plt.show()





#%%
input("Value를 저장하세요")
#%%

sds_wvlt = ap_wavelen
sds_wvlt.sort()
sds_val = np.array(ap_summed['aperture_sum']   )
sds_val = np.flip(sds_val)
#%%
"""
Data save
"""

Savepath=Path('/home/astro_02/AO2019-2/2019-10-24/Specdata')
#%%

table=Table()

wave_length = Column(data = sds_wvlt, name='Wavelength')
Value = Column(data = sds_val, name='Value')

table.add_column(wave_length, index = 0)
table.add_column(Value, index = 1)

table.write(Savepath/'NGC676_01.csv',format='ascii.csv',overwrite=True)


#%%
# Object만 바꿔서 저장용




#%%
Savepath=Path('/home/astro_02/AO2019-2/2019-10-24/Specdata')
DISPAXIS = 2 # 1 = line = python_axis_1 // 2 = column = python_axis_0
FONTSIZE = 12 # Change it on your computer if you wish.
rcParams.update({'font.size': FONTSIZE})
COMPIMAGE = ppdpath/'Comp-master.fits' # Change directory if needed!
LINE_FITTER = LevMarLSQFitter()

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

#for NGC2639
#ORDER_APTRACE = 9
#SIGMA_APTRACE = 9
#ITERS_APTRACE = 10
# The fitting is done by SIGMA_APTRACE-sigma ITERS_APTRACE-iters clipped on the
# residual of data
#%%
#for L in range(9):

for L in range(4):

#"""Change object here"""

    OBJIMAGE  = ppdpath/ f'WISEA-000{L+1}.fits'
  
    objhdu = fits.open(OBJIMAGE)

    objimage  = objhdu[0].data

    if DISPAXIS == 2:

        objimage = objimage.T
    elif DISPAXIS != 1:
        raise ValueError('DISPAXIS must be 1 or 2 (it is now {:d})'.format(DISPAXIS))
    
    EXPTIME = objhdu[0].header['EXPTIME']
    OBJNAME = objhdu[0].header['OBJECT']
    GAIN = objhdu[0].header['EGAIN']
    # Now python axis 0 (Y-direction) is the spatial axis 
    # and 1 (X-direciton) is the wavelength (dispersion) axis.
    N_SPATIAL, N_WAVELEN = np.shape(objimage)
    
    N_REID = N_SPATIAL//STEP_REID # No. of reidentification
    N_AP = N_WAVELEN//STEP_AP # No. of aperture finding
    
    # ``peak_local_max`` calculates the peak location using maximum filter:
    #   med1d_max = scipy.ndimage.maximum_filter(med1d, size=10, mode='constant')
    # I will use this to show peaks in a primitive manner.
    MINSEP_PK = 5   # minimum separation of peaks
    MINAMP_PK = 0.01 # fraction of minimum amplitude (wrt maximum) to regard as peak
    NMAX_PK = 50
    print("setting done!")
    
    
    # =============================================================================
    # apall (1): Plot a cut
    # =============================================================================
    lower_cut = N_WAVELEN//2 - NSUM_AP//2 
    upper_cut = N_WAVELEN//2 + NSUM_AP//2
    apall_1 = np.sum(objimage[:, lower_cut:upper_cut], axis=1)
    max_intens = np.max(apall_1)
    
    peak_pix = peak_local_max(apall_1, indices=True, num_peaks=10,
                              min_distance=int(N_SPATIAL/100),
                              threshold_abs=np.median(apall_1))
    
    x_apall = np.arange(0, len(apall_1))
    
    if L == 0:    
        fig = plt.figure()
        ax = fig.add_subplot(111)
        title_str = r'Peak ($A \geq {:.2f} A_\mathrm{{max}}$, $\Delta x \geq {:.0f}$)'
        
        ax.plot(x_apall, apall_1, lw=1)

    if L == 0:    
        for i in peak_pix:
            ax.plot((i, i), 
                    (apall_1[i]+0.01*max_intens, apall_1[i]+0.05*max_intens),
                    color='r', ls='-', lw=1)
            ax.annotate(i[0], (i, 0.9),
                        xycoords = ('data', 'axes fraction'),
                        fontsize='xx-small', rotation=70)
        ax.grid(ls=':')
        ax.set_xlabel('Pixel number')
        ax.set_ylabel('Pixel value')
        ax.set_xlim(0, len(apall_1))
        ax.set_title(title_str.format(np.median(apall_1), int(N_SPATIAL/100)))
        plt.show()
    
    
    # =============================================================================
    # apall(2): manually select sky, see how the fit works
    # =============================================================================
    ap_init = int(peak_pix[0])
    
    #for original sky region(HD207673,HD4628)
    #ap_sky = np.array([ap_init-40, ap_init-20, ap_init+20, ap_init+40])
    #for NGC 2639
    #ap_sky = np.array([ap_init-55, ap_init-35, ap_init+20, ap_init+40])
    #for WISEA
    ap_sky = np.array([ap_init-40, ap_init-20, ap_init+35, ap_init+55])
    
    # Regions to use as sky background. xl1 - 1, xu1, xl2 - 1, xu2. (0-indexing)
    #   Sky region should also move with aperture center!
    #   from ``ap_center - 50`` to ``ap_center - 40``, for example, should be used.
    # TODO: accept any even-number-sized ap_sky.
    
    # Interactive check
    x_sky = np.hstack( (np.arange(ap_sky[0], ap_sky[1]), 
                        np.arange(ap_sky[2], ap_sky[3])))
    sky_val = np.hstack( (apall_1[ap_sky[0]:ap_sky[1]], 
                          apall_1[ap_sky[2]:ap_sky[3]]))
    
    if L == 0:    
        fig = plt.figure(figsize=(8, 5))
        ax = fig.add_subplot(111)
        title_str = r'Skyfit: {:s} order {:d} ({:.1f}-sigma {:d}-iters)'
        ax.plot(x_apall, apall_1, lw=1)
        
    if FITTING_MODEL_APSKY.lower() == 'chebyshev':
        # TODO: maybe the clip is "3-sigma clip to residual and re-fit many times"?
        clip_mask = sigma_clip(sky_val, sigma=SIGMA_APSKY, iters=ITERS_APSKY).mask
        coeff_apsky, fitfull = chebfit(x_sky[~clip_mask], 
                                       sky_val[~clip_mask],
                                       deg=ORDER_APSKY,
                                       full=True)
        fitRMS = np.sqrt(fitfull[0][0]/n_found)
        n_sky = len(x_sky)
        n_rej = np.count_nonzero(clip_mask)
        sky_fit = chebval(x_apall, coeff_apsky) 
        
        if L == 0:

            ax.plot(x_apall, sky_fit, ls='--',
                    label='Sky Fit ({:d}/{:d} used)'.format(n_sky - n_rej, n_sky))
            ax.plot(x_sky[clip_mask], sky_val[clip_mask], marker='x', ls='', ms=10)
            [ax.axvline(i, lw=1, color='k', ls='--') for i in ap_sky]
            ax.legend()
            ax.set_title(title_str.format(FITTING_MODEL_APSKY, ORDER_APSKY,
                                          SIGMA_APSKY, ITERS_APSKY))
            ax.set_xlabel('Pixel number')
            ax.set_ylabel('Pixel value sum')
        
    else:
        raise ValueError('Function {:s} is not implemented.'.format(FITTING_MODEL_REID))
    if L == 0:
        ax.grid(ls=':')
        ax.set_xlabel('Pixel number')
        ax.set_ylabel('Pixel value')
        plt.show()
        
        
    # =============================================================================
    # apall (3): aperture trace
    # =============================================================================
    # within +- 100 pixels around the aperture, the wavelength does not change much
    # as can be seen from reidentify figure 
    # (in NHAO case, ~ 0.01% ~ 0.1 Angstrom order).
    # So it's safe to assume the wavelength is constant over around such region,
    # (spatial direction) and thus I will do sky fitting from this column,
    # without considering the wavelength change along a column.
    # Then aperture extraction will map the pixel to wavelength using aperture
    # trace solution.
    
    aptrace = []
    aptrace_fwhm = []
    #coeff_apsky = []
    #aptrace_apsum = []
    #aptrace_wavelen = []
    # TODO: This is quite slow as for loop used: improvement needed.
    # I guess the problem is sigma-clipping rather than fitting process..
    for i in range(N_AP - 1):
        lower_cut, upper_cut = i*STEP_AP, (i+1)*STEP_AP
        
        apall_i = np.sum(objimage[:, lower_cut:upper_cut], axis=1)
        sky_val = np.hstack( (apall_i[ap_sky[0]:ap_sky[1]], 
                              apall_i[ap_sky[2]:ap_sky[3]]))
    
        # Subtract fitted sky
        if FITTING_MODEL_APSKY.lower() == 'chebyshev':
            # TODO: maybe we can put smoothing function as IRAF APALL's b_naverage 
            clip_mask = sigma_clip(sky_val, sigma=SIGMA_APSKY, iters=ITERS_APSKY).mask
            coeff, fitfull = chebfit(x_sky[~clip_mask], 
                                     sky_val[~clip_mask],
                                     deg=ORDER_APSKY,
                                     full=True)
            apall_i -= chebval(x_apall, coeff)
    #        fitRMS = np.sqrt(fitfull[0][0]/n_found)
    #        n_sky = len(x_sky)
    #        n_rej = np.count_nonzero(clip_mask)
        
        else:
            raise ValueError('Function {:s} is not implemented.'.format(FITTING_MODEL_APSKY))
    
        #TODO: put something like "lost_factor" to multiply to FWHM_ID in the bounds.
        search_min = int(np.around(ap_init - 3*FWHM_AP))
        search_max = int(np.around(ap_init + 3*FWHM_AP))
        cropped = apall_i[search_min:search_max]
        x_cropped = np.arange(len(cropped))
        peak_pix = peak_local_max(cropped, 
                                  min_distance=FWHM_AP,
                                  indices=True,
                                  num_peaks=1)
        if len(peak_pix) == 0:
            aptrace.append(np.nan)
            continue
        peak_pix = peak_pix[0][0]
        
        #TODO: put something like "lost_factor" to multiply to FWHM_ID in the bounds.
        g_init = Gaussian1D(amplitude=cropped[peak_pix], 
                           mean=peak_pix, 
                           stddev=FWHM_AP * gaussian_fwhm_to_sigma,
                           bounds={'amplitude':(0, 2*cropped[peak_pix]) ,
                                   'mean':(peak_pix-3*FWHM_AP, peak_pix+3*FWHM_AP),
                                   'stddev':(0, FWHM_AP)})
        fitted = fitter(g_init, x_cropped, cropped)
        center_pix = fitted.mean.value + search_min
        std_pix = fitted.stddev.value
        aptrace_fwhm.append(fitted.fwhm)
        aptrace.append(center_pix)
    #    coeff_apsky.append(coeff)
    #    aptrace_apsum.append(apsum)
    #    apsum_lower = int(np.around(center_pix - apsum_sigma_lower * std_pix))
    #    apsum_upper = int(np.around(center_pix + apsum_sigma_upper * std_pix))
    #    apsum = np.sum(apall_i[apsum_lower:apsum_upper])
    
    aptrace = np.array(aptrace)
    aptrace_fwhm = np.array(aptrace_fwhm)
    #coeff_apsky = np.array(coeff_apsky)
    #aptrace_apsum = np.array(aptrace_apsum)
    
    
    
    # =============================================================================
    # apall(4): aperture trace fit
    # =============================================================================
    mask = np.isnan(aptrace)
    
    x_aptrace = np.arange(N_AP-1) * STEP_AP
    coeff_aptrace = chebfit(x_aptrace[~mask], aptrace[~mask], deg=ORDER_APTRACE)
    
    #aptrace_finite = aptrace[np.isfinite(aptrace)]
    #aptrace_finite
    
    resid_mask = sigma_clip(aptrace - chebval(x_aptrace, coeff_aptrace), 
                            sigma=SIGMA_APTRACE, iters=ITERS_APTRACE).mask
    
    x_aptrace_fin = x_aptrace[~resid_mask]
    aptrace_fin = aptrace[~resid_mask]
    coeff_aptrace_fin = chebfit(x_aptrace_fin, aptrace_fin, deg=ORDER_APTRACE)
    fit_aptrace_fin   = chebval(x_aptrace_fin, coeff_aptrace_fin)
    resid_aptrace_fin = aptrace_fin - fit_aptrace_fin
    del_aptrace = ~np.in1d(x_aptrace, x_aptrace_fin) # deleted points
    
    if L == 0:    
        fig = plt.figure(figsize=(10,8))
        gs = gridspec.GridSpec(3, 1)
        ax1 = plt.subplot(gs[0:2], sharex=ax2)
        ax2 = plt.subplot(gs[2])
        
        title_str = ('Aperture Trace Fit ({:s} order {:d})\n'
                    + 'Residuials {:.1f}-sigma, {:d}-iters clipped')
        plt.suptitle(title_str.format(FITTING_MODEL_APTRACE, ORDER_APTRACE,
                                      SIGMA_APTRACE, ITERS_APTRACE))
        ax1.plot(x_aptrace, aptrace, ls='', marker='+', ms=10)
        ax1.plot(x_aptrace_fin, fit_aptrace_fin, ls='--',
                 label="Aperture Trace ({:d}/{:d} used)".format(len(aptrace_fin), N_AP-1))
        ax1.plot(x_aptrace[del_aptrace], aptrace[del_aptrace], ls='', marker='x', ms=10)
        ax1.legend()
        ax2.plot(x_aptrace_fin, resid_aptrace_fin, ls='', marker='+')
        #ax2.plot(x_aptrace, aptrace - chebval(x_aptrace, coeff_aptrace_fin), 
        #         ls='', marker='+')
        ax2.axhline(+np.std(resid_aptrace_fin, ddof=1), ls=':', color='k')
        ax2.axhline(-np.std(resid_aptrace_fin, ddof=1), ls=':', color='k', 
                    label='residual std')
        
        ax1.set_ylabel('Found object position')
        ax2.set_ylabel('Residual (pixel)')
        ax2.set_xlabel('Dispersion axis (pixel)')
        ax1.grid(ls=':')
        ax2.grid(ls=':')
        ax2.set_ylim(-.5, .5)
        ax2.legend()
        plt.show()
    #plt.savefig('aptrace.png', bbox_inches='tight')
    
    
    # =============================================================================
    # apall(5): aperture sum
    # =============================================================================
#for original
#    apsum_sigma_lower = 5 # See below
#    apsum_sigma_upper = 5 
    #
#for WISEA,NGC2639
    apsum_sigma_lower = 3 # See below
    apsum_sigma_upper = 3 

    # lower and upper limits of aperture to set from the center in gauss-sigma unit.
    ap_fwhm = np.median(aptrace_fwhm[~resid_mask])
    ap_sigma = ap_fwhm * gaussian_fwhm_to_sigma
    apheight = (apsum_sigma_lower + apsum_sigma_upper)*ap_sigma
    
    x_ap = np.arange(N_WAVELEN)
    y_ap = chebval(x_ap, coeff_aptrace_fin)
    ap_wavelen = fit2D_REID(x_ap, y_ap)
    ap_sky_offset = ap_sky - ap_init
    
    ap_summed  = []
    data_skysub = []
    
    for i in range(N_WAVELEN):
        cut_i = objimage[:, i].copy()
        ap_sky_i = int(y_ap[i]) + ap_sky_offset
       
        x_sky = np.hstack( (np.arange(ap_sky_i[0], ap_sky_i[1]),
                            np.arange(ap_sky_i[2], ap_sky_i[3])))
        sky_val = np.hstack( (cut_i[ap_sky_i[0]:ap_sky_i[1]],
                              cut_i[ap_sky_i[2]:ap_sky_i[3]]))
        clip_mask = sigma_clip(sky_val, sigma=SIGMA_APSKY, maxiters=ITERS_APSKY).mask
       
        coeff = chebfit(x_sky[~clip_mask],
                        sky_val[~clip_mask],
                        deg=ORDER_APSKY)
    
        data_skysub.append(cut_i - chebval(np.arange(cut_i.shape[0]), coeff))
    
    data_skysub = np.array(data_skysub).T
    hdr = objhdu[0].header
    hdr.add_history(f"Sky subtracted using sky offset = {ap_sky_offset}, "
                    + f"{SIGMA_APSKY}-sigma {ITERS_APSKY}-iter clipping "
                    + f"to fit order {ORDER_APSKY} Chebyshev")
    _ = fits.PrimaryHDU(data=data_skysub, header=hdr)
    _.data = _.data.astype('float32')
    _.writeto(SKYSBPATH/(OBJIMAGE.stem+".skysub.fits"), overwrite=True)
    
    pos = np.array([x_ap, y_ap]).T
    aps = RectangularAperture(positions=pos, w=1, h=apheight, theta=0)
    phot = aperture_photometry(data_skysub, aps, method='subpixel', subpixels=30)
    ap_summed = phot['aperture_sum'] / EXPTIME

    if L == 0:
        
        fig, axs = plt.subplots(2, 1, figsize=(10, 6), sharex=False, sharey=False, gridspec_kw=None)
        
        axs[0].imshow(objimage, vmin=0, vmax=200, origin='lower')
        axs[1].imshow(data_skysub, vmin=0, vmax=200, origin='lower')
        axs[1].plot(pos[:,0], pos[:,1], 'r-')
        axs[0].set(title="Before sky subtraction")
        axs[1].set(title="Sky subtracted")
        
        for ax in [axs[0], axs[1]]:
            ax.plot(x_ap, y_ap+apsum_sigma_upper*ap_sigma, 'r-', lw=1)
            ax.plot(x_ap, y_ap-apsum_sigma_lower*ap_sigma, 'r-', lw=1)
        #    ax.set(ylim=(400, 700))
            for offset in ap_sky_offset:
                ax.plot(x_ap, y_ap+offset, 'y-', lw=1)
            
        plt.tight_layout()
        plt.show()   
        
        
        fig = plt.figure(figsize=(10,5))
        gs = gridspec.GridSpec(3,1)
        ax = plt.subplot(gs[0:3])
        
        #plt.setp(ax.get_xticklabels(), visible=False)
        title_str = ('{:s}\nobject = {:s}, EXPTIME = {:.1f} s')
        label_str = r'aperture width $\approx$ {:.1f} pix'
        plt.suptitle(title_str.format(OBJIMAGE.name,
                                      OBJNAME, EXPTIME))
        ax.plot(ap_wavelen, ap_summed, lw = 1,
                 label=label_str.format(apheight))
        ax.set_ylabel('Instrument Intensity\n(apsum/EXPTIME)')
        ax.set_xlabel(r'$\lambda$ ($\AA$)')
        
        ax.legend()
        plt.show()


        
    sds_wvlt = ap_wavelen
    sds_wvlt.sort()
    sds_val = np.array(ap_summed)
    sds_val = np.flip(sds_val)

    
    
    table=Table()
    
    wave_length = Column(data = sds_wvlt, name='Wavelength')
    Value = Column(data = sds_val, name='Value')
    
    table.add_column(wave_length, index = 0)
    table.add_column(Value, index = 1)
#"""Change object here"""
    
    table.write(Savepath/f'WISEA_raw_0{L+1}.csv',format='ascii.csv',overwrite=True)

