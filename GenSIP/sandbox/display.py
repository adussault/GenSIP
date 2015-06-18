"""
Contains functions for saving and displaying data used in histogram analysis
"""
import matplotlib.pyplot as plt
import numpy as np
import GenSIP.sandbox.datatools as dat
import GenSIP.sandbox.Histograms as hist
import os


def saveHist(Data,path,name='foil'):
    FIG=plt.figure()
    PLT = FIG.add_subplot(111)
    FIG.suptitle(name+" Histograms")
    
    plt.xlabel("Gray Level")
    plt.ylabel("Count")
    plotAllSects(Data,PLT)
    for reg in Data:
        if Data[reg]!='Null':
            Pt = hist.selectPtThresh(Data[reg])
            PLT.plot(Pt, dat.smoothed(Data[reg]['Histogram'])[Pt]+150,'cv')
    #Set Max and min of the y axis:
    assert os.path.isdir(path),"Not a directory."
    if not(os.path.exists(path+'/Histograms')):
        os.mkdir(path+'/Histograms')
    FIG.savefig(path+'/Histograms/'+name+".png")
    print "Histogram plot for "+name+" finished."
    plt.close()
   
def plotAllSects(Data,PLT, val='none'):
    labels = Data.keys()
    #labels = ['blk','pleat','darkMo','Mo','highEx','Plat']
    colors = ['black','purple','magenta','orange','green','cyan']
    for i in range(len(labels)):
        if (Data[labels[i]] != 'Null') and (type(Data[labels[i]]==dict)):
            print labels[i]
            plotFromDict(Data,labels[i],PLT,smooColor = colors[i])
            
def plotFromDict(Data, region, PLT, smooColor = 'orange'):
    '''
    This function recieves a dictionary of the format of an individual region 
    prodcued by MakeRegions and plots the histogram and histogram data of that region.
    Inputs:
        - Data - The data dictionary containing information on all regions and made
            by MakeRegions in the newmethod module.
        - region - a string giving the region. Must be identical to the key for
            the region in the Data dictionary.
    '''
    sectData = Data[region]
    PtThresh = hist.selectPtThresh(sectData)
    DirtThresh = hist.selectDirtThresh(sectData)
    x = np.arange(0,256,1)
    PLT.plot(x, dat.smoothed(sectData['Histogram']))
    PLT.plot(x, dat.nsmooth(sectData['Histogram'],i=6),smooColor)
    PLT.plot(sectData['Max'],dat.smoothed(sectData['Histogram'])[sectData['Max']],'ro')
    PLT.plot(sectData['Min'],dat.smoothed(sectData['Histogram'])[sectData['Min']],'ro')
    PLT.plot(sectData['Mean'],dat.smoothed(sectData['Histogram'])[sectData['Mean']],'go')
    PLT.plot(int(sectData['MoPeak']),dat.smoothed(sectData['Histogram'])[int(sectData['MoPeak'])],'y^')
    PLT.plot(PtThresh,dat.smoothed(sectData['Histogram'])[PtThresh]+dat.smoothed\
    (sectData['Histogram'])[int(sectData['MoPeak'])]/3,marker='v',color=smooColor,ms=10)
    PLT.plot(DirtThresh,dat.smoothed(sectData['Histogram'])[DirtThresh]+dat.smoothed\
    (sectData['Histogram'])[int(sectData['MoPeak'])]/2,marker='o',color=smooColor,ms=10)
    if sectData['Peaks'].size!=0:
        PLT.plot(sectData['Peaks'],dat.smoothed(sectData['Histogram'])[sectData['Peaks']],'g^')
    if sectData['Valleys'].size!=0:
        PLT.plot(sectData['Valleys'],dat.smoothed(sectData['Histogram'])[sectData['Valleys']],'bv')
    if sectData['NegInfl'].size!=0:   
        PLT.plot(sectData['NegInfl'],dat.smoothed(sectData['Histogram'])[sectData['NegInfl']],'rx')
    if sectData['PosInfl'].size!=0:
        PLT.plot(sectData['PosInfl'],dat.smoothed(sectData['Histogram'])[sectData['PosInfl']],'bx')
          
def RGBtoMoPtCrack(img):
    Mo = hist.threshMo(img).astype(np.uint8)
    Pt = hist.threshPt(img).astype(np.uint8)
    Crack = hist.threshCrack(img).astype(np.uint8)
    ret = np.zeros((img.shape[0],img.shape[1],3),dtype=np.uint8)
    ret[...,0]=Crack*255
    ret[...,1]=Pt*255
    ret[...,2]=Mo*255
    return ret
    
def RGBdisplay(imgs):
    
    ret = np.zeros((imgs[0].shape[0],imgs[0].shape[1],3),dtype=np.uint8)
    img2 = imgs[1].copy()
    img2[img2!=0]=1
    ret[...,0]=imgs[1]
    ret[...,1]=imgs[0]-img2*255
    ret[...,2]=imgs[0]-img2*255
    """
    for i in range(1,len(imgs)):
        imgs[i][imgs[i]!=0]=1
        imgs[i]=np.bitwise_not(imgs[i])
        try:
            ret[...,i]=ret[...,i]*imgs[i]
        except:
            pass
    """
    return ret