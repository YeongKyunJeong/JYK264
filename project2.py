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
for i in cptb['FILE']:
    cc = CCDData.read(i)
    images.append(cc)
mcomp = combine(images, method = 'median')
mcomp.header = cc.header
mcomp.header.add_history(f"{len(cptb)} image(s) median combined Comparision Ramp Image)")
mcomp = yfu.CCDData_astype(mcomp,dtype = 'float32')
mcomp.write(ppdpath/'Comp-master.fits',overwrite=True)

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
print(Path(OBJIMAGE).exists())
