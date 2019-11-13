import os
import glob
from astropy.io import fits
from astropy.table import Table, Column
import numpy as np
from pathlib import Path
from matplotlib import pyplot as plt
import matplotlib.colors as colors

import ysfitsutilpy as yfu

import astropy.units as u
from astropy.nddata import CCDData
from ccdproc import Combiner, combine
from ccdproc import subtract_bias, subtract_dark, flat_correct
from ccdproc import trim_image


#%%
# 아래에 fits파일이 있는 directory 위치를 적을 것.
# 원본 파일 경로와 preprocess를 할 작업 파일 경로 설정 및 directory 변경

fitspath = Path('/home/astro_02/AO2019-2/2019-10-24/BU')
newfitspath = Path('/home/astro_02/AO2019-2/2019-10-24')
os.chdir(fitspath)
allfitsname = glob.glob('*.fits')
allfitspath = list(fitspath.glob('*.fit'))

print(len(allfitspath),"Items are searched")


#%%
# Image Trimming 여부와 범위 설정, Default

df = input("Default? Y/N(Otherwise) [880:1100,1,2048]")
if df != "Y":

    print("이미지를 자르시겠습니까?")
    cutting_decision = input("Y/N, 그 외의 모든 입력은 N으로 처리")

    if cutting_decision == "Y":
        x1= int(input("x시작"))
        x2= int(input("x끝"))
        y1= int(input("y시작"))
        y2= int(input("y끝"))
        print(f"[{x1}:{x2}, {y1}:{y2}]")
    # 
    if not ( (x1 < x2) & (y1 < y2)):
        print("범위 입력이 잘못되었습니다. Image trimming을 하지 않습니다.")
        cutting_decision = "N"
        
    else:
        x1 = 1
        x2 = 2048
        y1 = 1
        y2 = 2048
        
else:
    x1 = 880
    x2 = 1100
    y1 = 1
    y2 = 2048
    print(f"[{x1}:{x2},{y1}:{y2}]")


#%%
# 이미지 trimming과 작업 파일 복사, 복사 후 chmod 777 *.fits을 터미널에 입력할 것

"""header test"""
hdr = fits.getheader(allfitspath[0])
print("hdr을 입력해서 header를 읽어보시오")

# Table Variables Setting
cards = ['DATE-OBS', 'NAXIS1', 'NAXIS2',
         'XPIXSZ','YPIXSZ','XBINNING', 'YBINNING',
         'EXPTIME', 'EGAIN',
         'OBJECT']

dtypes = ['U24', int, int,
          float,float, int, int,
          float, float,
          'U16']

df = input("현재 작업 파일을 그대로 사용하겠습니까? Y/N, 그 외의 모든 입력은 N으로 처리")

if df == "Y":
    for filename in allfitsname:
        ccd = trim_image(CCDData.read(filename,unit=u.adu),fits_section=f"[{x1}:{x2}, {y1}:{y2}]")
        #16bit를 32비트로 만듦
        ccd = yfu.CCDData_astype(ccd,dtype = 'float32')
        filename += 's'
        ccd.write(newfitspath/filename,overwrite=True)    
        #fits로 만들어 저장
    
    
#%%
# Trimming된 새 Data들로 Path와 directory를 바꿈
    
os.chdir(newfitspath)
allfitsname = glob.glob('*.fits')
allfitspath = list(newfitspath.glob('*.fits'))
allfitspath.sort()
allfitsname.sort()


#%%
# 과제 1번의 분류를 위해 fdic으로 분류했으나 아래에서는 table만 사용하므로 하지 않아도 무관

table=Table(names=cards, dtype=dtypes)
fnames=[]
obnames=[]
fdic = dict(bias=[], dark={}, flat=[], comp=[], objt=[])

for filename in allfitsname:
    fnames.append(filename)
    hdr = fits.getheader(filename)
    row = []
    for card in cards:
        row.append(hdr[card])
    table.add_row(row)
    
    if hdr['OBJECT'].lower().startswith('calibration'):
        try:
            fdic['dark'][float(hdr['EXPTIME'])]
        except KeyError:
            fdic['dark'][float(hdr['EXPTIME'])] = []
        fdic['dark'][float(hdr['EXPTIME'])].append(filename)        
        
    elif hdr['OBJECT'].lower().startswith('bias'):
        fdic['bias'].append(filename)
    
    elif hdr['OBJECT'].lower().startswith('flat'):
        fdic['flat'].append(filename)
    elif hdr['OBJECT'].lower().startswith('comp'):
        fdic['comp'].append(filename)
    else:
        fdic['objt'].append(filename)
        obnames.append(hdr["OBJECT"])

# Table에 이름 Column 추가  
fnames = Column(data=fnames, name='FILE')
table.add_column(fnames, index = 0)
table.sort('FILE')

# 작업 파일 Directory 아래에 Headertable이라는 Directory를 만들어 Table을 따로 저장
tablepath=Path(newfitspath/'Headertable')
if not tablepath.exists():
    tablepath.mkdir()
else:
    print("Table 저장 폴더가 이미 있어서 새로 만들지 않음")
    
table.write(tablepath/'Headersumm.csv',format='ascii.csv',overwrite=True)


#%%
# Preprocess Image Table Setting

biastab = table[(table['OBJECT'] == 'Bias')]

# Exposure Time List를 만들어 해당 index에 대응되는 Exposure time을 부를 때 사용
# Dark는 dictionary를 사용해서 Exposure Time별로 Table 생성
darkdic={}
exptlist=[]
for i in table['EXPTIME']:
    if ((not i in exptlist ) & (i != 0)):
        exptlist.append(int(i))
exptlist.sort()
for i in range(len(exptlist)):
    dark = table[((table['OBJECT']=='Calibration') & (table['EXPTIME'] == exptlist[i]))]
    darkdic[exptlist[i]]=dark

flattab = table[(table['OBJECT'] == 'Flat')]


#%%
# Preprocess Image Name Setting, Dark는 Dictionary로 Exposure Time별로 따로 입력

bias_fname = 'bias.fits'

dark_fdic = {}
for i in range(len(exptlist)):
    dark_fdic[exptlist[i]] = f'dark{exptlist[i]}s.fits'
    
flat_fname = 'flat.fits'


#%%
# 작업 파일 Directory 아래에 Preprocess이라는 Directory를 만들어 Table을 따로 저장

prepath=Path(newfitspath/'Preprocess')
if not prepath.exists():
    prepath.mkdir()
else:
    print("Prepath 저장 폴더가 이미 있어서 새로 만들지 않음")

if os.path.exists(prepath/bias_fname):
    mbias = CCDData.read(prepath/bias_fname, unit = u.adu)
    print("이전에 만든 bias 사용")
else:   
    images=[]
    # Bias Image를 불러 Median Combine
    for i in range(len(biastab)):
        cc = CCDData.read(biastab[i]['FILE'], unit = u.adu)
        images.append(cc)
    mbias = combine(images, method = 'median')
    
    # 첫 번째 Raw Flat의 Header를 가져와 Master Flat에 넣고 History 추가, 32bit로 바꿔 저장
    mbias.header = cc.header
    mbias.header.add_history(f"{len(biastab)} image(s) median combined bias frames")
    mbias = yfu.CCDData_astype(mbias,dtype = 'float32')
    mbias.write(prepath/bias_fname,overwrite=True)
    print(type(mbias))

fig, axs = plt.subplots(1,1)
im = yfu.zimshow(axs,mbias)
axs.set_xlabel("Bias", fontsize=12)

mbias = CCDData.read(prepath/bias_fname, unit = u.adu)


#%%
# Master Dark이 Data를 Dictionary로 저장

mdark_fdic={}
for k in range(len(exptlist)):
    if os.path.exists(prepath/dark_fdic[exptlist[k]]):
        mdark_fdic[exptlist[k]] = CCDData.read(prepath/dark_fdic[exptlist[k]]
                                    ,unit = u.adu)
        print(f"이전에 만든 {exptlist[k]}초 Dark 사용")

    else:   
        images=[]
        for i in range(len(darkdic[exptlist[k]])):
            cc = CCDData.read(darkdic[exptlist[k]][i]['FILE'], unit= u.adu)
            images.append(cc)
            print(darkdic[exptlist[k]][i]['FILE'])
        print(len(images), "장의 Image를 Combine")
        mdark_fdic[exptlist[k]] = combine(images,method = 'median')
        print(f"{exptlist[k]}초 Master Dark",sep='\n')
        
        # Bias Subtraction        
        mdark_fdic[exptlist[k]] = subtract_bias(mdark_fdic[exptlist[k]],mbias)
        
        # Header 추가
        mdark_fdic[exptlist[k]].header = cc.header
        mdark_fdic[exptlist[k]].header.add_history(f"{len(darkdic[exptlist[k]])} image(s) median combined dark frames with bias sub")
        mdark_fdic[exptlist[k]] = yfu.CCDData_astype(mdark_fdic[exptlist[k]],dtype = 'float32')
        mdark_fdic[exptlist[k]].write(prepath/dark_fdic[exptlist[k]],overwrite=True) 
        
fig, axs = plt.subplots(1,len(exptlist))
for i in range(len(exptlist)):
    im = yfu.zimshow(axs[i],mdark_fdic[exptlist[i]])
    axs[i].set_xlabel(f"{exptlist[i]} s", fontsize=12)

mdark_fdic={}
for k in range(len(exptlist)):
    if os.path.exists(prepath/dark_fdic[exptlist[k]]):
        mdark_fdic[exptlist[k]] = CCDData.read(prepath/dark_fdic[exptlist[k]]
                                    ,unit = u.adu)

    
#%%

if os.path.exists(prepath/flat_fname):
    mflat = CCDData.read(prepath/flat_fname, unit = u.adu)
    print("이전에 만든 bias 사용")    
    
else:
    images=[]
    for i in range(len(flattab)):
        cc = CCDData.read(flattab[i]['FILE'], unit = u.adu)
        cc = subtract_bias(cc,mbias)
        cc = subtract_dark(cc,mdark_fdic[int(flattab[i]['EXPTIME'])]
                            ,flattab[i]['EXPTIME']*u.second,flattab[i]['EXPTIME']*u.second)
        # 각 Image의 평균으로 Normalize
        hdr = cc.header
        cc = CCDData(cc/np.mean(cc),unit = u.adu)
        cc.header = hdr
        images.append(cc)
    mflat = combine(images,method = 'median')
    
    mflat.header = cc.header
    mflat.header.add_history(f"{len(flattab)} image(s) median combined flat frames with b&d sub")
    mflat = yfu.CCDData_astype(mflat,dtype = 'float32')
    mflat.write(prepath/flat_fname,overwrite=True)

fig, axs = plt.subplots(1,1)
im = yfu.zimshow(axs,mflat)
axs.set_xlabel("Flat", fontsize=12)

mflat = CCDData.read(prepath/flat_fname, unit = u.adu)


#%%
# Object Image Table Setting

ob_dic={}
ob_list=[]
# Cal도 Flat도 Bias도 아니면 Object, Dictionay로 object별로 저장
for i in table['OBJECT']:
    if ((not i in ob_list ) & 
        (i !=  'Calibration') & 
        (i != 'Flat') & 
        (i != 'Bias') ):
        ob_list.append(i)
ob_list.sort()
for i in range(len(ob_list)):
    ob = table[((table['OBJECT'] == ob_list[i]))]
    ob_dic[ob_list[i]] = ob


#%%
# Image를 Data와 Header를 따로 저장
    
ppdpath = Path(newfitspath/'PreprocessedImages')
if not ppdpath.exists():
    ppdpath.mkdir()

for i in ob_dic:
    print(i)
    for j in range(len(ob_dic[i])):
        c = CCDData(fits.getdata(ob_dic[i][j]['FILE']), unit = u.adu)
        hdr = fits.getheader(ob_dic[i][j]['FILE'])
        #Exposure Time에 맞는 Dark
        mdark = mdark_fdic[hdr['EXPTIME']]
        c = subtract_bias(c,mbias)
        c = subtract_dark(c,mdark
                            ,hdr['EXPTIME']*u.second,hdr['EXPTIME']*u.second)
        c = flat_correct(c,mflat,min_value=0.01)
        
        c.header = hdr
        c.header.add_history("B & D subtracted and F corrected")
        c = yfu.CCDData_astype(c,dtype = 'float32')
        c.write(ppdpath/ob_dic[i][j]['FILE'],overwrite=True) 
        
        
#%%
        
fig, axs = plt.subplots(2,len(ob_list),figsize=(12,12))
for i in range(len(ob_list)):
    data1 = CCDData.read(f"{ob_list[i]}-0001.fits")
    data2 = CCDData.read(ppdpath/f"{ob_list[i]}-0001.fits")
    im1 = yfu.zimshow(axs[0][i],data1)
    im2 = yfu.zimshow(axs[1][i],data2)
    axs[0][i].set_xlabel(f"{ob_list[i]}\nraw", fontsize=8)
    axs[1][i].set_xlabel(f"{ob_list[i]}", fontsize=8)














