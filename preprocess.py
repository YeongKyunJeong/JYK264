import os
import glob
from astropy.io import fits
from astropy.table import Table, Column
import numpy as np
from pathlib import Path
from matplotlib import pyplot as plt
import matplotlib.colors as colors

import astropy.units as u
from astropy.nddata import CCDData
from ccdproc import Combiner

#%%
# 아래에 fits파일이 있는 directory 위치를 적을 것.

#집용
#namepath = 'C:/Users/jjbyk/AO2/project/2019-10-24'
#fitspath = Path('C:/Users/jjbyk/AO2/project/2019-10-24')

#학과전산실용
namepath = '/home/astro_02/AO2019-2/2019-10-24'
fitspath = Path('/home/astro_02/AO2019-2/2019-10-24')


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


flattab = table[(table['OBJECT'] == 'Flat')]

""" 30초가 맞는지 10초가 맞는지 상의할 것"""






"""Preprocess Image fits 만들기 >>> 이미지 합치는 법을 알아와서 아래 else문 고칠 것"""

bias_fname = 'bias.fits'

dark_fdic = {}
for i in range(len(exptlist)):
    dark_fdic[exptlist[i]] = f'dark{exptlist[i]}s.fits'
    
flat_fname = 'flat.fits'






#%%
prepath=Path(fitspath/'Preprocessed')
if not prepath.exists():
    prepath.mkdir()
else:
    print("Prepath 저장 폴더가 이미 있어서 새로 만들지 않음")



if os.path.exists(bias_fname):
    mbias = fits.getdata(bias_fname)

else:   
    images=[]
    for i in range(len(biastab)):
        cc = CCDData(fits.getdata(biastab[i]['FILE']), unit = u.adu)
        images.append(cc)
        
    cc = Combiner(images)
    mbias = cc.median_combine()
    cc = fits.getheader(biastab[0]['FILE'])
    
    mbias.header = cc
    mbias.header.add_history(f"{len(biastab)} image(s) median combined bias frames")
    mbias.write(prepath/bias_fname,overwrite=True)




mdark_fdic={}
for i in range(len(exptlist)):
    if os.path.exists(dark_fdic[exptlist[i]]):
        mdark_fdic[exptlist[i]] = fits.getdata(dark_fdic[exptlist[i]])

else:
    for j in exptlist:   

        images=[]
        for i in range(len(darkdic[j])):
            cc = CCDData(fits.getdata(darkdic[j][i]['FILE']),unit= u.adu)
            images.append(cc)
        cc = Combiner(images)
        mdark_fdic[j] = cc.median_combine()
        cc = fits.getheader(darkdic[j][0]['FILE'])
        
        mdark_fdic[j].header = cc
        mdark_fdic[j].header.add_history(f"{len(darkdic[j])} image(s) median combined dark frames")
        mdark_fdic[j].write(prepath/dark_fdic[j],overwrite=True)    


#%%

if os.path.exists(flat_fname):
    mflat = fits.getdata(flat_fname) 

else:

    images=[]
    for i in range(len(flattab)):
        cc = CCDData(fits.getdata(flattab[i]['FILE']), unit = u.adu)
        images.append(cc)

    cc = Combiner(images)
    mflat = cc.median_combine()
    cc = fits.getheader(flattab[0]['FILE'])

    mflat.header = cc
    mflat.header.add_history(f"{len(flattab)} image(s) median combined flat frames")
    mflat.write(prepath/flat_fname,overwrite=True)

"""
else:
    mflat = preproc.make_master_flat(flattab, mbias=mbias, sigma=3, iters=5,
                                     min_value=5000, 
                                     output = flat_fname)
"""
#%%
mflat.header
#%%

"""Object Image Table 나누기"""


#%%
test1 = CCDData(np.ones((4,4)),unit=u.adu)
test2 = CCDData(np.zeros((4,4)),unit=u.adu)
test3 = CCDData(np.ones((4,4)),unit=u.adu)
c = Combiner([test1, test2, test3])
#%%
ccdavg = c.average_combine()
ccdmed = c.median_combine()
print(ccdavg,ccdmed,sep='\n')
print(type(test1))
#%%
fig = plt.figure(figsize=(10,5))
ax1 = fig.add_subplot(1,2,1)
ax2 = fig.add_subplot(1,2,2)
im1 = ax1.imshow(test1)
im2 = ax2.imshow(test1)
#%%
test = []
for i in range(len(biastab)):
    cc = CCDData(fits.getdata(biastab[i]['FILE']),unit= u.adu)
    test.append(cc)
print(len(test))
#test1 = CCDData(fits.getdata('Bias-0001.fit'),unit= u.adu)
#c = Combiner([test1])
#print(test1.shape,type(test1))
#%%
cc = Combiner(test)

#%%

ccdavg = cc.average_combine()
ccdmed = cc.median_combine()
print(ccdavg,ccdmed,sep='\n')


#%%
for i in range(500,510):

    testnum=[]
    for j in range(len(biastab)):
        testnum.append(float(fits.getdata(f'Bias-000{j+1}.fit')[i][i]))
    print(np.mean(testnum),np.median(testnum))
#%%
for i in range(500,510):
    print(ccdavg[i][i],ccdmed[i][i])
#실행용이므로 남겨둘 것















