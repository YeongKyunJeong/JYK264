
import os
import glob
from astropy.io import fits
from astropy.table import Table, Column
import numpy as np
from pathlib import Path
from matplotlib import pyplot as plt
import matplotlib.colors as colors

#%%
namepath = 'C:/Users/jjbyk/AO2/project/2019-10-24'
fitspath = Path('C:/Users/jjbyk/AO2/project/2019-10-24')
os.chdir(namepath)
#os.chdir("..") directory 변경
#allfits = fitspath.glob('*.fit')
# recursive=True로 설정하고 '**'를 사용하면 모든 하위 디렉토리까지 탐색한다.

"""Image 이름, 경로"""

allfitsname = glob.glob('*.fit')
allfitspath = list(fitspath.glob('*.fit'))

print(len(allfitsname),type(allfitsname[0]),"\n",len(allfitspath),type(allfitspath[0]))

#%%
"""
Header 읽기, Table 작성
"""
#%%
"""header test"""
hdr = fits.getheader(allfitsname[0])
"""hdr을 입력해서 header를 읽어보시오"""




"""Table 변수 정리"""
cards = ['DATE-OBS', 'NAXIS1', 'NAXIS2',
         'XPIXSZ','YPIXSZ','XBINNING', 'YBINNING',
         'EXPTIME', 'EGAIN',
         'OBJECT']

dtypes = ['U24', int, int,
          float,float, int, int,
          float, float,
          'U16']




table=Table(names=cards, dtype=dtypes)
fnames=[]

for filename in allfitsname:
    os.chmod(filename,777)
    fnames.append(filename)
    hdr = fits.getheader(filename)
    row = []
    for card in cards:
        row.append(hdr[card])
    table.add_row(row)    

# table에 이름 column 추가
fnames = Column(data=fnames, name='FILE')
table.add_column(fnames, index = 0)
table.sort('FILE')
# table 폴더를 따로 만들어 저장함
tablepath=Path(fitspath/'Headertable')
if not tablepath.exists():
    tablepath.mkdir()
else:
    print("Table 저장 폴더가 이미 있어서 새로 만들지 않음")
    
table.write(tablepath/'Headersumm.csv',format='ascii.csv',overwrite=True)





"""Preprocess용 Table 나누기"""
biastab = table[(table['OBJECT'] == 'Bias')]
# exptime 가짓 수
darkdic={}
exptlist=[]
for i in table['EXPTIME']:
    if ((not i in exptlist ) & (i != 0)):
        exptlist.append(int(i))
exptlist.sort()

#print(exptlist)

for i in range(len(exptlist)):
    dark = table[((table['OBJECT']=='Calibration') & (table['EXPTIME'] == exptlist[i]))]
    darkdic[exptlist[i]]=dark

#flat의 exptime이 하나일 경우
#flattab = table[(table['OBJECT'] == 'Flat')]

#이번 경우(expt = 10)만 사용
flattab = table[((table['OBJECT'] == 'Flat') & (table['EXPTIME'] == 10))]

#30초가 맞는지 10초가 맞는지, 둘 다 쓰는지 몰라 만든 예비
#flattab10 = table[((table['OBJECT'] == 'Flat') & (table['EXPTIME'] == 10))]
#flattab30 = table[((table['OBJECT'] == 'Flat') & (table['EXPTIME'] == 30))]

""" 30초가 맞는지 10초가 맞는지 상의할 것"""






"""Preprocess Image fits 만들기 >>> 이미지 합치는 법을 알아와서 아래 else문 고칠 것"""

bias_fname = 'bias.fits'

dark_fdic = {}
for i in range(len(exptlist)):
    dark_fdic[exptlist[i]] = f'dark{exptlist[i]}s.fits'
    
flat_fname = 'flat.fits'







if os.path.exists(bias_fname):
    mbias = fits.getdata(bias_fname)

"""
else:
    mbias = preproc.make_master_bias(biastab, sigma=3, iters=5, min_value=0,
                                     output = bias_fname)
"""
 
mdark_fdic={}
for i in range(len(exptlist)):
    if os.path.exists(dark_fdic[exptlist[i]]):
        mdark_fdic[exptlist[i]] = fits.getdata(dark_fdic[exptlist[i]])

"""        
    else:
         mdark_fdic[exptlist[i]] = preproc.make_master_dark(darkdic[exptlist[i]],
                                     mbias=mbias, sigma=3, iters=5,
                                     min_value=0, 
                                     output = dark_fname)
"""        
    

if os.path.exists(flat_fname):
    mflat = fits.getdata(flat_fname) 

"""
else:
    mflat = preproc.make_master_flat(flattab, mbias=mbias, sigma=3, iters=5,
                                     min_value=5000, 
                                     output = flat_fname)
"""
#%%

"""Object Image Table 나누기"""


#%%
print(darkdic[exptlist[0]])
print(dark_fdic[exptlist[0]])
print(mdark_fdic[exptlist[0]])

#실행용이므로 남겨둘 것















