"""
histogram_tools
Contains functions that identify features in the histogram of an image or region
in an image and determine the platinum or dirt threshold values. 

"""   

import numpy as np
import GenSIP.histomethod.datatools as dat

###################################################################################

###################################################################################

def selectDirtThresh(sectData, verbose=True):
    """selects the dirt threshold value given the data dictionary of a region
    in the dictionary produced by NewRegThresh"""
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
        if verbose: print "MoPeak: "+str(MoPeak)
        Pot_Dirt = np.arange(0,MoPeak,1)
        Pot_Dirt = Pot_Dirt[Histogram[Pot_Dirt]<=Histogram[MoPeak]/20]
        for i in range(1,21):
            if Pot_Dirt[Histogram[Pot_Dirt]<=Histogram[MoPeak]/(20/i)].size!=0:
                Dirt_Thresh = Pot_Dirt.max()
                break
        else:
            if verbose: print "Dirt thresh set to zero."
            Dirt_Thresh = 0
    #while n == 0:
    #    if Potential_Infl.max():
    #        pass
    return Dirt_Thresh
    
###################################################################################

###################################################################################

def selectPtThresh(sectData, verbose=True):
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
    
    if verbose: print "MoPeak: "+str(sectData['MoPeak'])
    
    if PtPeaks.size!=0:
        # Filter out any peaks that may be close to the MoPeak
        PtPeaks = PtPeaks[sectData['Histogram'][PtPeaks]<0.9*sectData['Histogram'][sectData['MoPeak']]]
    if PtPeaks.size!=0:
        if verbose: print "Potential Pt peaks found."
        PtPeak = PtPeaks[sectData['Histogram'][PtPeaks].argmax()]
        PtVals = PtVals[PtVals<PtPeak] # Remove any 
        if PtVals.size != 0:
            PtThresh = np.max(PtVals)
            if verbose: print "Pt Threshold value set to a valley between the Molybdenum peak and Pt peak."
        else:
            PtThresh = PtPeak
            if verbose: print "Pt Threshold value set to the Platinum peak."
    else:
        if verbose: print "No potential Pt peaks."
        if PtInfls.size != 0:
            PtThresh = np.median(PtInfls)
            if verbose: 
                print "Pt Threshold value set to the median inflection point between the MoPeak and PtPeak."
        elif PtVals.size != 0:
            PtThresh = np.max(PtVals)
            # The smoothing algorithms actually move the location of the valley
            # Find the actually valley by looking for the minimum in the Histogram
            # near the PtThresh in the smoothed Histogram.
            PtThresh = sectData['Histogram'][PtThresh-4:PtThresh+5].argmin()+PtThresh-4
            if verbose: print "Pt Threshold set to Valley between MoPeak and the max value."
        else:
            PtThresh = int(sectData['Max'])
            if verbose: print "Pt Threshold set to maximum pixel value in the image."
    if PtThresh<100:
        if verbose: print "Pt Thresh set to 100 to prevent it from being set too low."
        PtThresh = 100
    if verbose: print "Pt Thresh: " + str(PtThresh)
    return int(PtThresh)

###################################################################################

###################################################################################

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


        