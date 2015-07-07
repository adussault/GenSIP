"""
_________________
HISTOGRAMS MODULE\______________________________________________________________

"""   

import numpy as np
import GenSIP.histomethod.datatools as dat


def selectDirtThresh(sectData):
    Histogram = sectData['Histogram']
    MoPeak = sectData['MoPeak']
    Pot_Infl = sectData['PosInfl'][sectData['PosInfl']<MoPeak]
    n = 0
    # Filter the Potential inflection points so that they are at maximum one tenth
    # the number of counts as the molybdenum peak. 
    Pot_Infl = Pot_Infl[Histogram[Pot_Infl]<=Histogram[MoPeak]/10]
    # Make sure the Inflection point is 
    if (Pot_Infl.size!=0) and (Histogram[Pot_Infl.max()] >= Histogram[MoPeak]/25):
        Dirt_Thresh = Pot_Infl.max()
    else:
        print "MoPeak: "+str(MoPeak)
        Pot_Dirt = np.arange(0,MoPeak,1)
        Pot_Dirt = Pot_Dirt[Histogram[Pot_Dirt]<=Histogram[MoPeak]/20]
        for i in range(1,21):
            if Pot_Dirt[Histogram[Pot_Dirt]<=Histogram[MoPeak]/(20/i)].size!=0:
                Dirt_Thresh = Pot_Dirt.max()
                break
        else:
            print "Dirt thresh set to zero."
            Dirt_Thresh = 0
    #while n == 0:
    #    if Potential_Infl.max():
    #        pass
    return Dirt_Thresh
    
def selectPtThresh(sectData):
    """
    Takes the data dictionary of a specific region and uses the Histogram to determine
    the optimum threshold value for the platinum. 
    """
    # All peaks with higher grayscale vals than MoPeak
    PtPeaks = sectData['Peaks'][sectData['Peaks']>sectData['MoPeak']] 
    
    # All valleys with higher grayscale vals than MoPeak
    PtVals = sectData['Valleys'][sectData['Valleys']>sectData['MoPeak']] 
    # All negative Inflection points with higher grayscale vals than MoPeak
    PtInfls = sectData['NegInfl'][sectData['NegInfl']>sectData['MoPeak']] 
        
    if PtPeaks.size!=0:
        # Filter out any peaks that may be close to the MoPeak
        PtPeaks = PtPeaks[sectData['Histogram'][PtPeaks]<0.9*sectData['Histogram'][sectData['MoPeak']]]
    if PtPeaks.size!=0:
        PtPeak = PtPeaks[sectData['Histogram'][PtPeaks].argmax()]
        PtVals = PtVals[PtVals<PtPeak] # Remove any 
        if PtVals.size != 0:
            PtThresh = np.max(PtVals)
        else:
            PtThresh = PtPeak
    else:
        if PtInfls.size != 0:
            PtThresh = np.median(PtInfls)
        elif PtVals.size != 0:
            PtThresh = np.max(PtVals)
            # The smoothing algorithms actually move the location of the valley
            # Find the actually valley by looking for the minimum in the Histogram
            # near the PtThresh in the smoothed Histogram.
            PtThresh = sectData['Histogram'][PtThresh-4:PtThresh+5].argmin()+PtThresh-4
        else:
            PtThresh = int(sectData['Max'])
    if PtThresh<100:
        PtThresh = 100
    return int(PtThresh)
    
def IDPeaks(img):
    #PtMax = np.max(img)
    #DirtCrackMin = np.min(img)
    MoX, histo = findMoPeakByImg(img)
    #x = np.arange(256)
    PksX,PksY = dat.getMaxima(histo)
    ValX,ValY = dat.getMinima(histo)
    PksX,PksY = np.asarray(PksX),np.asarray(PksY)
    ValX,ValY = np.asarray(ValX),np.asarray(ValY)
    
    dif = [abs(i-MoX) for i in PksX]
    dif = np.asarray(dif)
    # ID Molybdenum Peak
    MoPkX = PksX[dif.argmin()]
    #MoPkY = histo[MoX]
    
    # Platinum Peak assumed to be the largest peak with a higher graylevel value than the Mo peak
    
    PtX,PtY = PksX[PksX>MoPkX],PksY[PksX>MoPkX]
    PtPkX = PksX[PksY.argmax()]
    
    # Find Valley between the Mo and Pt Peaks
    MoPtX,MoPtY = ValX[(MoPkX<ValX)&(ValX<PtPkX)],ValY[(MoPkX<ValX)&(ValX<PtPkX)]
    #MoPtValX = MoPtX[MoPtY.argmin()]
    
    
    
    
"""
_______________________
THRESHOLDING OPERATIONS\________________________________________________________

"""      
    
        
def threshMo(img):
    peakVal, histo = findMoPeakByImg(img)
    MIN,MAX = dat.atFWHM(histo,peakVal)
    thresh = dat.easyThresh(img,MIN,MAX,Binary=True)
    return thresh

def threshPt(img):
    return (img==255)
    
def threshCrack(img):
    return (img==0)

def findMoPeakByImg(img, returnHist = False):
    """
    Takes the histogram of an image and determines the likely maximum of the 
    Molybdenum peak. Returns the pixel value (0 to 255) of the maximum of the 
    molybdenum peak.
    # The Molybdenum peak is usually the biggest and widest peak, but not the tallest.
    # Therefore, I first take the histogram of the image with 32 sections of 8 pixel-values
    # each, and then I determine which of those 32 sections has the highest count value.
    # Then I look for the maximum value within that section in a histogram divided into 256 
    parts, and that is the peak value. 
    """
    histo = np.histogram(img,256,(0,255))
    histo8 = np.histogram(img,32,(0,255))
    counts = histo[0]
    #Smoothed = dat.smoothed(counts)
    
    #widths = np.ones(Smoothed.shape)
    
    M = histo8[0].argmax()*8
    firstPeak = counts[M:M+8].argmax()+M
    firstPeak = int(firstPeak)
    if returnHist:
        return firstPeak,counts
    else:
        return firstPeak

'''
def findMoPeakByHistogram(histo):
    """
    Takes the histogram of an image and determines the likely maximum of the 
    Molybdenum peak. Returns the pixel value (0 to 255) of the maximum of the 
    molybdenum peak.
    """
    #Darkest = min(img)
    #Lightest = max(img)
    histo = np.histogram(img,256,(0,255))
    histo8 = np.histogram(img,32,(0,255))
    counts = histo[0]
    #Smoothed = dat.smoothed(counts)
    
    #widths = np.ones(Smoothed.shape)
    
    M = histo8[0].argmax()*8
    firstPeak = counts[M:M+8].argmax()+M
    firstPeak = int(firstPeak)
    return firstPeak,counts
'''

        