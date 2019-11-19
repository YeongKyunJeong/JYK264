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
from ccdproc import Combiner
from ccdproc import subtract_bias, subtract_dark, flat_correct
from ccdproc import trim_image

#%%
# 아래에 fits파일이 있는 directory 위치를 적을 것.

#집용
#namepath = 'C:/Users/jjbyk/AO2/project/2019-10-24'
#fitspath = Path('C:/Users/jjbyk/AO2/project/2019-10-24')

#학과전산실용
fitspath = Path('/home/astro_02/AO2019-2/2019-10-24/BU')
newfitspath = Path('/home/astro_02/AO2019-2/2019-10-24')

os.chdir(fitspath)
#os.chdir("..") directory 변경
#allfits = fitspath.glob('*.fit')
# recursive=True로 설정하고 '**'를 사용하면 모든 하위 디렉토리까지 탐색한다.

"""Image 이름, 경로"""
allfitsname = glob.glob('*.fits')
allfitspath = list(fitspath.glob('*.fit'))


print(len(allfitspath),"Items are searched")

#%%
"""
Header 읽기, Table 작성
"""

#%%
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
    
    if not ( (x1 < x2) & (y1 < y2)):
        
        print("범위 입력이 잘못되었습니다. Image trim을 하지 않습니다.")
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


"""header test"""
hdr = fits.getheader(allfitspath[0])
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

"""Image trimming and coppying"""


for filename in allfitsname:
    ccd = trim_image(CCDData.read(filename,unit=u.adu),fits_section=f"[{x1}:{x2}, {y1}:{y2}]")
    #16bit를 32비트로 만듦
    ccd = yfu.CCDData_astype(ccd,dtype = 'float32')
    filename += 's'
    ccd.write(newfitspath/filename,overwrite=True)    
    #fits로 만들어 저장
    
    
#%%
os.chdir(newfitspath)
allfitsname = glob.glob('*.fits')
allfitspath = list(newfitspath.glob('*.fits'))
allfitspath.sort()
#Trim된 새 Data들로 Path를 바꿈
#%%

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


# table에 이름 column 추가  
fnames = Column(data=fnames, name='FILE')
table.add_column(fnames, index = 0)
table.sort('FILE')
# table 폴더를 따로 만들어 저장함
tablepath=Path(newfitspath/'Headertable')
if not tablepath.exists():
    tablepath.mkdir()
else:
    print("Table 저장 폴더가 이미 있어서 새로 만들지 않음")
    
table.write(tablepath/'Headersumm.csv',format='ascii.csv',overwrite=True)


#%%



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
#%%






"""Preprocess Image fits 만들기 >>> 이미지 합치는 법을 알아와서 아래 else문 고칠 것"""

bias_fname = 'bias.fits'

dark_fdic = {}
for i in range(len(exptlist)):
    dark_fdic[exptlist[i]] = f'dark{exptlist[i]}s.fits'
    
flat_fname = 'flat.fits'






#%%
prepath=Path(newfitspath/'Preprocessed')

if not prepath.exists():
    prepath.mkdir()
else:
    print("Prepath 저장 폴더가 이미 있어서 새로 만들지 않음")



if os.path.exists(prepath/bias_fname):
    mbias = fits.getdata(prepath/bias_fname)
    print("이전에 만든 bias 사용")
else:   
    images=[]
    for i in range(len(biastab)):
        cc = CCDData(fits.getdata(biastab[i]['FILE']), unit = u.adu)
        
        cc = yfu.CCDData_astype(cc,dtype = 'float32')
        images.append(cc)
        
    cc = Combiner(images)

    mbias = cc.median_combine()
    cc = fits.getheader(biastab[0]['FILE'])
    
    mbias.header = cc
    mbias.header.add_history(f"{len(biastab)} image(s) median combined bias frames")
    mbias = yfu.CCDData_astype(mbias,dtype = 'float32')
    mbias.write(prepath/bias_fname,overwrite=True)
#%%
mdark_fdic={}
for k in range(len(exptlist)):
    if os.path.exists(prepath/dark_fdic[exptlist[k]]):
        mdark_fdic[exptlist[k]] = fits.getdata(prepath/dark_fdic[exptlist[k]])
        print(f"이전에 만든 {exptlist[k]}초 Dark 사용")

    else:   
        images=[]
        for i in range(len(darkdic[exptlist[k]])):
            cc = CCDData(fits.getdata(darkdic[exptlist[k]][i]['FILE']),unit= u.adu)
            images.append(cc)
            print(darkdic[exptlist[k]][i]['FILE'])

        print(len(images))
        cc = Combiner(images)
        mdark_fdic[exptlist[k]] = cc.median_combine()
        cc = fits.getheader(darkdic[exptlist[k]][0]['FILE'])
        print(f"{exptlist[k]}초 Master Dark",sep='\n')
        
        mdark_fdic[exptlist[k]] = subtract_bias(mdark_fdic[exptlist[k]],mbias)
        #bias 빼기
        mdark_fdic[exptlist[k]].header = cc
        mdark_fdic[exptlist[k]].header.add_history(f"{len(darkdic[exptlist[k]])} image(s) median combined dark frames with bias sub")
        mdark_fdic[exptlist[k]] = yfu.CCDData_astype(mdark_fdic[exptlist[k]],dtype = 'float32')
        mdark_fdic[exptlist[k]].write(prepath/dark_fdic[exptlist[k]],overwrite=True) 


     
#%%

if os.path.exists(prepath/flat_fname):
    mflat = fits.getdata(prepath/flat_fname)
    print("이전에 만든 bias 사용")
# flat에 해당하는 expt가 없을 경우는 나중에 추가
else:

    images=[]
    for i in range(len(flattab)):
        cc = CCDData(fits.getdata(flattab[i]['FILE']), unit = u.adu)
#        print(np.mean(cc))
        cc = subtract_bias(cc,mbias)
        #bias 빼기
        cc = subtract_dark(cc,mdark_fdic[int(flattab[i]['EXPTIME'])]
                            ,flattab[i]['EXPTIME']*u.second,flattab[i]['EXPTIME']*u.second)
        #dark 빼기
        cc = CCDData(cc/np.mean(cc),unit = u.adu)
        #Normalization

        images.append(cc)

    cc = Combiner(images)
    mflat = cc.median_combine()
    cc = fits.getheader(flattab[0]['FILE'])

    mflat.header = cc
    mflat.header.add_history(f"{len(flattab)} image(s) median combined flat frames with b&d sub")
    mflat = yfu.CCDData_astype(mflat,dtype = 'float32')
    mflat.write(prepath/flat_fname,overwrite=True)

#%%

"""Object Image Table 나누기"""
ob_dic={}
ob_list=[]
for i in table['OBJECT']:
    if ((not i in ob_list ) & 
        (i !=  'Calibration') & 
        (i != 'Flat') & 
        (i != 'Bias') ):
        ob_list.append(i)
ob_list.sort()

#print(exptlist)

for i in range(len(ob_list)):
    ob = table[((table['OBJECT']==ob_list[i]))]
    ob_dic[ob_list[i]]=ob

#%%
i
for i in ob_dic:
    print(i)
    for j in range(len(ob_dic[i])):
        c = CCDData(fits.getdata(ob_dic[i][j]['FILE']), unit = u.adu)
        hdr = fits.getheader(ob_dic[i][j]['FILE'])
        mdark = mdark_fdic[hdr['EXPTIME']]
        
        
        c = subtract_bias(c,mbias)
        c = subtract_dark(c,mdark
                            ,hdr['EXPTIME']*u.second,hdr['EXPTIME']*u.second)
        c = flat_correct(c,mflat,min_value=0.01)
        
        c.header = hdr
        c.header.add_history("B & D subtracted and F corrected")
        c = yfu.CCDData_astype(c,dtype = 'float32')
        c.write(prepath/ob_dic[i][j]['FILE'],overwrite=True) 

#%%
"""
#%%
c = CCDData(fits.getdata(ob_dic['Comp'][0]['FILE']), unit = u.adu)
hdr = fits.getheader(ob_dic['Comp'][0]['FILE'])
#%%
prepath/mdark_fdic[ str( hdr['EXPTIME'])]
#%%
#print(mdark_fdic[10])

print(np.mean(mbias),np.mean(mdark))

#%%
hdr


print(ob_list, ob_dic,sep='\n')
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
fig, axs = plt.subplots(1,1)
im = yfu.zimshow(axs,mbias)
axs.set_xlabel("Bias", fontsize=12)


#%%
fig, axs = plt.subplots(1,len(exptlist))
for i in range(len(exptlist)):
    im = yfu.zimshow(axs[i],mdark_fdic[exptlist[i]])
    axs[i].set_xlabel(f"{exptlist[i]} s", fontsize=12)

#%%
#fig, axs = plt.subplots(1, len(exptlist))
#st = yfu.give_stats(mdark_fdic[exptlist[-1]].data)
#vmin = st["zmin"]
#vmax = st["zmax"]
#for i in range(len(exptlist)):
#    im = axs[i].imshow(mdark_fdic[exptlist[i]], origin='lower', vmin=vmin, vmax=vmax)
#    axs[i].set_xlabel
    


#%%

for i in range(500,510):
    print(ccdavg[i][i],ccdmed[i][i])
#실행용이므로 남겨둘 것

"""













