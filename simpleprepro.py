import ysfitsutilpy as yfu
from pathlib import Path

top = Path("2019-10-24")
ngc2639paths = top.glob("NGC2639*.fits")
preprocdir = top/"Preprocess"

mbiaspath = preprocdir/"bias.fits"
mflatpath = preprocdir/"flat.fits"
outdir = Path("test")
#%%
for fpath in ngc2639paths:
    ccd_xxx = yfu.load_ccd(fpath)
    exptime = ccd_xxx.header["EXPTIME"]
    ccd_bdx = yfu.bdf_process(ccd_xxx,
                              mbiaspath=mbiaspath,
                              mdarkpath=preprocdir/f"dark{exptime:.0f}s.fits",
                              output=outdir/f"{fpath.stem}_bdx.fits")
    ccd_bdf = yfu.bdf_process(ccd_xxx,
                              mbiaspath=mbiaspath,
                              mdarkpath=preprocdir/f"dark{exptime:.0f}s.fits",
                              mflatpath=mflatpath,
                              output=outdir/f"{fpath.stem}_bdf.fits")
    
    