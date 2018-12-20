#!/usr/bin/python
import matplotlib
matplotlib.use('agg')
import csv
import pylab
from pyteomics import fasta, parser,mzml,mgf
from pyteomics import auxiliary as aux
import pandas as pd
import glob
import itertools
import seaborn as sns
import numpy as np
from statsmodels.nonparametric.smoothers_lowess import lowess
import os
from pyteomics import pylab_aux
import matplotlib.pyplot as plt
import argparse

def read_data(name_mzml, name_mzml_calibration, name_psm,delim,colname):
    print ("reading data")
    injtime_ms1=[]
    injtime_ms2=[]
    starttime_ms1=[]
    starttime_ms2=[]
    indexms1=[]
    charge_ms2=[]
    mz_ms2=[]
    angle_x=[]
    angle_y=[]
    prec_int=[]
    psm=set()
    x=[]
    y=[]
    a=0
    for sc in mzml.read(name_mzml):
        if sc['ms level']==1:
            injtime_ms1.append(float(sc['scanList']['scan'][0]['ion injection time']))
            starttime_ms1.append(float(sc['scanList']['scan'][0]['scan start time']))
            indexms1.append(int(sc['id'].split(" ")[2].split('=')[1]))
        if sc['ms level']==2:
            it_ms2=float(sc['scanList']['scan'][0]['ion injection time'])
            st_ms2=float(sc['scanList']['scan'][0]['scan start time'])
            try:
                charge=int(sc["precursorList"]['precursor'][0]['selectedIonList']['selectedIon'][0]['charge state'])
            except (KeyError, ValueError):
                charge=0
            mz=sc["precursorList"]['precursor'][0]['selectedIonList']['selectedIon'][0]['selected ion m/z']
            try:
                prec_intensity=sc["precursorList"]['precursor'][0]['selectedIonList']['selectedIon'][0]['peak intensity']
            except (KeyError, ValueError):
                a+=1
            if sc['intensity array'].size:
                intensity=sc['intensity array'].sum(dtype=float)/sc['intensity array'].size
            else:
                intensity=0
            prec_int.append(prec_intensity)
            injtime_ms2.append(it_ms2)
            angle_y.append(intensity)
            angle_x.append(len(sc['intensity array']))
            starttime_ms2.append(st_ms2)
            mz_ms2.append(mz)
            charge_ms2.append(charge)

    coef=10**(len(str(int(np.mean(angle_y))))-1)

    #read_for_angle_calibration
    if name_psm is not None:
        print ("reading reference files for angle score calculation")
        with open(name_psm) as fIn:
            reader = csv.DictReader(fIn, delimiter=delim)
            next(reader)
            for row in reader:
                psm.add(int(row[colname].split('.')[1]))
        for sc in mzml.read(name_mzml_calibration):
            if sc['ms level']==2 and len(sc["m/z array"])!=0 and int(sc['id'].split(' ')[2].split('=')[1])+form in psm:
                x1=len(sc["m/z array"])
                x.append(x1)
                y1=float(sum(sc['intensity array']))/(len(sc['intensity array']))
                y.append(y1)
        coef=10**(len(str(int(np.mean(y))))-1)
        y=np.array(y)/coef
    print ("reading is finished, calculating metrics")
    return (injtime_ms1, injtime_ms2, starttime_ms1, starttime_ms2, indexms1, charge_ms2, mz_ms2, angle_x, angle_y, prec_int, psm, x, y,coef)

def ms1_ms2(starttime_ms1,starttime_ms2,ax):
    width = 0.5
    ms1=len(starttime_ms1)
    ms2=len(starttime_ms2)
    rects1 = ax.bar([0],ms1, width,alpha=1,color=color[0])
    plt.text(0, ms1/2, ms1, ha='center', fontsize=20)
    rects1 = ax.bar([1],ms2, width,alpha=1,color=color[1])
    plt.text(1, ms2/2, ms2, ha='center', fontsize=20)
    ax.set_xticks(np.arange(2))
    xTickMarks = ['MS1','MS2']
    ax.set_xticks(np.arange(2))
    xtickNames = ax.set_xticklabels(xTickMarks)
    plt.setp(xtickNames, rotation=0, fontsize=15)
    ax.set_title('MS1/MS2')

def aqtime(starttime_ms1,starttime_ms2):
    starttime_all=sorted([(x,1) for x in starttime_ms1]+[(x,2) for x in starttime_ms2])
    aqtime_ms1=[]
    aqtime_ms2=[]
    for i,k in enumerate(starttime_all[:-1]):
        if k[1]==1 and any( [(k[1]-starttime_all[i+1][1])==-1,(k[1]-starttime_all[i+1][1])==0] ):
            aqtime_ms1.append(starttime_all[i+1][0]-k[0])
        if k[1]==2 and any( [(k[1]-starttime_all[i+1][1])==1,(k[1]-starttime_all[i+1][1])==0] ):
            aqtime_ms2.append(starttime_all[i+1][0]-k[0])
    a=pylab.hist(np.array(aqtime_ms1)*60 ,histtype='step', lw=2,density=True, label='MS1, sum=%.2f sec'%sum(aqtime_ms1),color=color[0])
    b=pylab.hist(np.array(aqtime_ms2)*60 ,histtype='step', lw=2, density=True, label='MS2, sum=%.2f sec'%sum(aqtime_ms2),color=color[1])
    pylab.legend(loc=1,fontsize=12)
    pylab.xlabel("AT, sec",fontsize=15)
    pylab.ylabel("# spectra, normalized", fontsize=15)
    pylab.title('Acquisition time')

def it_ms1(starttime_ms1,injtime_ms1):
    pylab.scatter(starttime_ms1, injtime_ms1,s=10, alpha=0.7,color=color[0])
    pylab.legend(markerscale=2)
    pylab.xlabel("Start Time, min",fontsize=15)
    pylab.ylabel("Injection Time,ms", fontsize=15)
    pylab.title('MS1')

def inten_prec(starttime_ms2,start,finish,prec_int):
    ind=np.logical_and(np.array(starttime_ms2)>start, np.array(starttime_ms2)<finish)
    prec=np.log10(np.array(prec_int))[ind]
    b=np.percentile(prec,99.9)
    pylab.hist(prec,bins=np.linspace(0,max(prec),100), color=color[0],alpha=0.5, lw=1,edgecolor='k' )
    pylab.xlabel("Log10(intensity)",fontsize=15)
    pylab.xlim(np.percentile(prec,0.1),b)
    pylab.ylabel("Scans #",fontsize=15)
    pylab.title('Intensity precursor ions')

def it_ms2(starttime_ms2,start,finish,injtime_ms2):
    ind=np.logical_and(np.array(starttime_ms2)>start, np.array(starttime_ms2)<finish)
    pylab.hist(np.array(injtime_ms2)[ind],bins=np.linspace(0,max(injtime_ms2),100), color=color[1],alpha=0.5, lw=1,edgecolor='k')
    pylab.xlabel("Injection Time,ms",fontsize=15)
    pylab.ylabel("Scans #",fontsize=15)
    pylab.title('MS2_IT')

def realtop(starttime_ms1,indexms1,):
    pylab.scatter(np.array(starttime_ms1)[:-1], np.ediff1d(indexms1)-1, alpha=0.7, s=15, label='TopN',color=color[1])
    fit=lowess(np.ediff1d(indexms1)-1,np.array(starttime_ms1)[:-1], frac=0.05, it=0)
    pylab.plot(fit[:,0],fit[:,1], "r-" ,label='Average TopN')
    pylab.ylim(0, max(np.ediff1d(indexms1))+max(np.ediff1d(indexms1))/10)
    pylab.xlim(0, max(starttime_ms1)+max(starttime_ms1)/15)
    pylab.legend(loc=1,markerscale=2,fontsize=15)
    pylab.xlabel('RT, min',fontsize=15)
    pylab.ylabel('TopN',fontsize=15)
    pylab.title('Real_TopN',fontsize=15)

def charge(maxcharge,charge_ms2,starttime_ms2,mz_ms2):
    for i in range(0,maxcharge+1):
        charge_ms2=np.array(charge_ms2)
        starttime_ms2=np.array(starttime_ms2)
        mz_ms2=np.array(mz_ms2)
        msk=[k==i for k in charge_ms2]
        pylab.scatter(np.array(starttime_ms2)[msk], np.array(mz_ms2)[msk], s=20, alpha=0.5,color=color[i-1], label='charge_%s+(%i scans)'%(i,len(np.array(mz_ms2)[msk])))
        pylab.title('Precursor ions',fontsize=15)
        pylab.legend(markerscale=1.5,fontsize=15)
        pylab.xlabel('RT,min',fontsize=15)
        pylab.ylabel('m/z',fontsize=15)

def angle_calculation(x,y,angle_x,angle_y,coef):
    x_mod=[]
    y_mod=[]
    per_1_x=np.percentile(x,2)
    per_99_x=np.percentile(x,98)
    per_1_y=np.percentile(y,2)
    per_99_y=np.percentile(y,98)
    for i in zip(x,y):
        if per_1_x<i[0]<per_99_x and per_1_y<i[1]<per_99_y:
            x_mod.append(i[0])
            y_mod.append(i[1])
    for i in np.arange(0.001,1,0.001):
        dif=[]
        for k in zip(x_mod,y_mod):
            dif.append(i*(k[0]-per_1_x)+per_1_y-k[1])
        under_line=sum(x > 0 for x in dif)
        if float(under_line)/len(x_mod)>0.01:
            under=i
            break
        else:
            continue
    for j in np.arange(1,0.001,-0.01):
        dif=[]
        for k in zip(x_mod,y_mod):
            dif.append(j*(k[0]-per_1_x)+per_1_y-k[1])
        above_line=sum(x < 0 for x in dif)
        if float(above_line)/len(x_mod)>0.01:
            above=j
            break
        else:
            continue
    a=0
    angle_y_mod=np.array(angle_y)/coef
    for x1,y1 in zip(angle_x,angle_y_mod):
        if y1-(x1-per_1_x)*above<per_1_y and y1-(x1-per_1_x)*under>per_1_y:
            a=a+1
    per=100*float(a)/len(angle_x)
    return (under,above,per_1_x,per_1_y,per,angle_y_mod,coef)

def angle (under,above,per_1_x,per_1_y,per,angle_x,angle_y_mod,coef):
    a=np.arange(per_1_x,1000,100.)
    b=under*(a-per_1_x)+per_1_y
    pylab.plot(a,b,color='black')
    c=above*(a-per_1_x)+per_1_y
    pylab.plot(a,c,color='black')
    pylab.scatter(angle_x,angle_y_mod, alpha=0.5, color=color[1],label='score=%.2f'%per)
    pylab.xlim(0,np.percentile(angle_x,99.5))
    pylab.ylim(0,np.percentile(angle_y_mod,99.5))
    pylab.xlabel('#peaks',fontsize=15)
    pylab.legend(fontsize=15, loc=1)
    pylab.ylabel('avg intensity, 10^%i'%np.log10(coef),fontsize=15)

def inten_number_peaks_ms1(angle_x,angle_y):
    coef=10**(len(str(int(np.mean(angle_y))))-1)
    angle_y_mod=np.array(angle_y)/coef
    pylab.scatter(angle_x,angle_y_mod, alpha=0.5, color=color[1])
    pylab.xlim(0,np.percentile(angle_x,99.5))
    pylab.ylim(0,np.percentile(angle_y_mod,99.5))
    pylab.xlabel('#peaks',fontsize=15)
    pylab.ylabel('avg intensity, 10^%i'%np.log10(coef),fontsize=15)
    pylab.title('MS/MS',fontsize=15)

#Arguments reading

parser = argparse.ArgumentParser()

parser.add_argument('input', help='mzML file with path')
parser.add_argument('-o', '--output', help='path to save result, by default save in the folder of the input file')
parser.add_argument('-refPSM', nargs='?', help='csv file with psm identifications for angle score calculating')
parser.add_argument('-refmzML', nargs='?', help='mzML file for angle score calculating')
parser.add_argument('-d', nargs='?', help='delimiter in csv file with psm identifications for angle score calculating; tab by default',default='\t')
parser.add_argument('-cn', nargs='?', help='column name with spectra names in csv file with psm identifications for angle score calculating; "spectrum" by default',default='spectrum')
parser.add_argument('-start', nargs='?', help='delay time before sample actually comes to mass spec; using for precursor intensity and injection time (MS/MS) calculation; 0 by default',default=0)
parser.add_argument('-stop', nargs='?', help='time of wash starting; using for precursor intensity and injection time (MS/MS) calculation. By default maximum analysis time')
parser.add_argument('-charge', nargs='?', help='max charge of precursor ions; 4 by default',default=4)
parser.add_argument('-f', nargs='?', help='format of input file for identification process, 0 for mzML and 1 for mgf; 1 by default',default=4)
args = parser.parse_args()

if os.path.exists(args.input):
    name_mzml=args.input
else:
    print("Could not find the input file %s" % args.input)
    sys.exit(1)

if args.output is None:
      output=os.path.split(name_mzml)[0]
else:
      output=args.output

if args.refPSM is None or args.refmzML is None:
    print ("There is no files for angle score calculation, continue without it")

name=os.path.split(name_mzml)[1].split('.')[0]
form=int(args.f)
name_psm=args.refPSM
name_mzml_calibr=args.refmzML
start=int(args.start)
maxcharge=int(args.charge)
delim=str(args.d)
colname=str(args.cn)

if (args.refPSM is not None) and (os.path.split(name_mzml_calibr)[1].split('.')[0] not in name_psm):
     print ("Warning! File names for angle score calibration don't match")

injtime_ms1, injtime_ms2, starttime_ms1, starttime_ms2, indexms1, charge_ms2, mz_ms2, angle_x, angle_y, prec_int, psm, x, y, coef= read_data(name_mzml,name_mzml_calibr,name_psm,delim,colname)

if args.stop is None:
    finish=max(starttime_ms1)
else:
    finish=int(args.stop)


#for pretty pictures
plt.style.use('seaborn-whitegrid')
plt.rcParams['ytick.labelsize'] = 15
plt.rcParams['xtick.labelsize'] = 15
plt.rcParams['axes.titlesize'] = 15
plt.rcParams['axes.titlesize']=15
color=["#b84c7d","#4d8ac9","#4bc490","#7f63b8",'b','g','#edd630']

#build pictures

pylab.figure(figsize=(15,40))
ax1 = plt.subplot2grid((6, 2), (0, 0))
ms1_ms2(starttime_ms1,starttime_ms2,ax1)
ax2 = plt.subplot2grid((6, 2), (0, 1))
aqtime(starttime_ms1,starttime_ms2)
ax3 = plt.subplot2grid((6, 2), (1, 0), colspan=2)
it_ms1(starttime_ms1,injtime_ms1)
ax4 = plt.subplot2grid((6, 2), (2, 0), colspan=2)
inten_prec(starttime_ms2,start,finish,prec_int)
ax5 = plt.subplot2grid((6, 2), (3, 0), colspan=2)
it_ms2(starttime_ms2,start,finish,injtime_ms2)
ax6 = plt.subplot2grid((6, 2), (4, 0), colspan=2)
charge(maxcharge,charge_ms2,starttime_ms2,mz_ms2)
ax7 = plt.subplot2grid((6, 2), (5, 0))
realtop(starttime_ms1,indexms1)
if name_psm is not None:
    under,above,per_1_x,per_1_y,per,angle_y_mod,coef=angle_calculation(x,y,angle_x,angle_y,coef)
    ax8 = plt.subplot2grid((6, 2), (5, 1))
    angle (under,above,per_1_x,per_1_y,per,angle_x,angle_y_mod,coef)
else:
    ax8 = plt.subplot2grid((6, 2), (5, 1))
    inten_number_peaks_ms1(angle_x,angle_y)
print ("saving results")

pylab.savefig('%s/viQC_results_%s.png'%(output,name))

print ('enjoy your QC')

