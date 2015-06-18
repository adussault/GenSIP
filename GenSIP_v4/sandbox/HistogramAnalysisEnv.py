# This is a script for setting up tools for analyzing histograms
import GenSIP_v4.sandbox.foldertools as f
import GenSIP_v4.sandbox.newmethod as nm
import GenSIP_v4.sandbox.Histograms as hist
import GenSIP_v4.sandbox.display as dis

import GenSIP_v4.functions as fun
import GenSIP_v4.bigscans.images as im

import os
import numpy as np
import random
import matplotlib.pyplot as plt

foilPath = "InputPicts/FoilScans/40360_2/sub_imgs_Q1"
maskPath = "InputPicts/FoilScans/40360_2/sub_imgs_Q1_mask"
names = []
Max={}
Min={}
MoPeaks = {}
Histos = {}
FWHMs = {}
Means = {}
decentSubs = []

x = np.arange(0,256,1)

preposts = {}
posters = {}
Data = {}
IMAGES = f.DoForFolder(foilPath,fun.loadImg, filter=f.FILonlySubimages,flag=0)
masks = f.DoForFolder(maskPath, fun.loadImg, filter=f.FILonlySubimages,flag=0)

for i in IMAGES:
    names.append(i)

names.sort()

for i in names:
    prepost=nm.PosterPreProc(IMAGES[i])
    posters[i] = im.bigPosterfy(prepost)

for i in names:
    Data[i] = nm.MakeRegions(IMAGES[i],posters[i], Mask=masks[i])

# select out the standards for special work
stds = f.FILonlySubimages(os.listdir("standards"))

# for plotting
Colours = {}
rand_color = lambda: random.randint(0,255)
rand_hexcolor = lambda: '#%02X%02X%02X' % (rand_color(),rand_color(),rand_color())
for i in names:
     colour = rand_hexcolor()
     Colours[i] = colour
FIG = plt.figure()
PLT = FIG.add_subplot(111)

for i in stds:
    for j in Data[i]:
        if Data[i][j] != 'Null':
            dis.plotFromDict(Data[i],j,PLT,smooColor = Colours[i])
'''
for i in names:
    Max[i] = np.max(IMAGES[i])

for i in names:
    Min[i] = np.min(IMAGES[i])

for i in names:
    MoPeaks[i],Histos[i] = findMo(IMAGES[i])

for i in names:
    FWHMs[i]=atFWHM(Histos[i],MoPeaks[i])
'''   
'''
AllPeaks = {}
AllValleys = {}

for i in names:
    smoo = smoothed(Histos[i])
    AllPeaks[i],Zy = getMaxima(smoo)
    AllValleys[i],Zy = getMinima(smoo)

for i in names:
    if 90<MoPeaks[i]:
        decentSubs.append(i)

selections=[]

for i in range(6):
    selections.append(decentSubs[random.randint(0,len(decentSubs)-1)])

    




OGselects = ['sub_004_002',
 'sub_013_010',
 'sub_010_011',
 'sub_010_015',
 'sub_004_007',
 'sub_006_006']

FIG = plt.figure()
plt1 = FIG.add_subplot(111)
    
for name in selections:
    smoo = smoothed(Histos[name])
    Zx,Zy = getMaxima(smoo, smoonum=4)
    cursmoo = np.gradient(np.gradient(smoo))*5
    plt1.plot(x,cursmoo,Colours[name])
    plt1.plot(MoPeaks[name],smoo[MoPeaks[name]],'r^')
    plt1.plot(Zx,Zy,'g^',linestyle='None')
    plt1.plot(AllValleys[name],smoo[AllValleys[name]],'bv')
    plt1.plot(x,smoo,Colours[name])



    
nameByIndex=np.ndarray((16,16))
'''
"""
PLOTTING HISTOGRAMS:
    
plot(x,Histos[decentsubs[6]],figure=FIG)
plot(Peaks[decentSubs[6]],Histos[decentSubs[6]][Peaks[decentSubs[6]]],'g^')
plot(FWHMs[decentSubs[6]][0],Histos[decentSubs[6]][FWHMs[decentSubs[6]][0]],'ro')
plot(FWHMs[decentSubs[6]][1],Histos[decentSubs[6]][FWHMs[decentSubs[6]][1]],'ro')

EXTRA CODE:
    
FINDING % Pt
for i in names:
    if Data[i]['Mo'] != 'Null':
        PtThresh = selectPtThresh(Data[i]['Mo'])
        print i+': '+str(PtThresh)
        print '   MoPeak: ' + str(Data[i]['Mo']['MoPeak'])
        print '   Max: ' + str(Data[i]['Mo']['Max'])
        print '% Pt: ' + str(int(100*float(np.sum(Data[i]['Mo']['Histogram'][Data[i]['Mo']['Histogram']>PtThresh]))/float(np.sum(Data[i]['Mo']['Histogram']))))
    else:
        print i + ' Null'
WRITING RGB PT IMAGES:
    for i in names:
    Pt = ThreshImage(IMAGES[i]).astype(np.uint8)
    ret = np.zeros((img.shape[0],img.shape[1],3),dtype=np.uint8)
    ret[...,1]=Pt
    ret[...,2]=IMAGES[i].astype(np.uint8)
    cv2.imwrite('Output/Output_40360_2_Q1/NewPtMaps/'+i+'.png',ret,[cv2.cv.CV_IMWRITE_PNG_COMPRESSION,6])
    
for i in names:
    Pt = ThreshImage(IMAGES[i]).astype(np.uint8)
    ret = np.zeros((IMAGES[i].shape[0],IMAGES[i].shape[1],3),dtype=np.uint8)
    ret[...,0]=Posters[i].astype(np.uint8)/2
    ret[...,1]=Pt
    ret[...,2]=IMAGES[i].astype(np.uint8)
    cv2.imwrite('Output/Output_40360_2_Q1/NewPtMaps/'+i+'.png',ret,[cv2.cv.CV_IMWRITE_PNG_COMPRESSION,6])
    saveHist(Data[i],i)


"""