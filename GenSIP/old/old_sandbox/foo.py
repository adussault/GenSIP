################################################################################
###     FOO is a temporary file for writing and testing new functions        ###
################################################################################

import matplotlib.pyplot as plt
import cv2
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.backends.backend_tkagg as tkagg
from matplotlib.figure import Figure
from scipy import misc
import GenSIP.functions as fun
import os
#from GenSIP.sandbox.Histograms import floodBySeed


"""
_____________________________
CREATING A MASK FOR BIG FOILS\__________________________________________________

"""

BWImgs = "InputPicts/test_imgs/particle_count"
PlatImgs = "InputPicts/test_imgs/plat_count/"
DirtImgs = "InputPicts/test_imgs/dirt_count/"

OGimgs = [m for m in os.listdir(DirtImgs) if m.endswith('.png')]
for img in OGimgs:
    #os.rename(BWImgs+'/'+img,BWImgs+'/'+img[0:28]+'.png')
    os.rename(PlatImgs+img,PlatImgs+img[0:28]+'.png')
    os.rename(DirtImgs+img,DirtImgs+img[0:28]+'.png')
'''
for img in OGimgs:
    platImg = fun.loadImg(BWImgs+'/'+img).astype(np.uint8)
    dirtImg = platImg.copy()
    platImg[platImg==platImg.max()]=230
    platImg[platImg!=platImg.max()]=185
    dirtImg[dirtImg==dirtImg.max()]=40
    dirtImg[dirtImg!=dirtImg.max()]=185
    cv2.imwrite(PlatImgs+img,platImg,[cv2.cv.CV_IMWRITE_PNG_COMPRESSION,6])
    cv2.imwrite(DirtImgs+img,dirtImg,[cv2.cv.CV_IMWRITE_PNG_COMPRESSION,6])
'''

    


'''
def bigMaskEdges(image,res, maxFeatureSize=2000, Bkgrdthreshold = 95):
    """
    Creates a mask for masking off the edges of a large foil scan.
        Inputs:
         - image - image to be masked
        Key-word Arguments:
         - maxFeatureSize - largest feature diameter not considered to be the 
            edge of the foil. Default set to 10000 microns. 
         - threshold - approximate maximum value of the background
    """
    # want to reduce image to an image with a height of 1 mm
    rszheight = int(1000/np.sqrt(res)) 
    print "Resize height: " + str(rszheight)
    rszfactor = float(rszheight)/float(image.shape[0])
    print "Resize factor: " + str(rszfactor)
    resized = misc.imresize(image, rszfactor, interp='bilinear')
    print "Resize shape: " + str(resized.shape)
    scaledMaxFeat = int(rszfactor*maxFeatureSize/np.sqrt(res))
    print "Kernel width: " + str(scaledMaxFeat)
    threshed = resized.astype(np.uint8)
    threshed[threshed<=Bkgrdthreshold] = 0
    threshed[threshed>Bkgrdthreshold] = 1
    
    kernel = np.ones((scaledMaxFeat,scaledMaxFeat))
    removeDirt = cv2.morphologyEx(threshed, cv2.MORPH_CLOSE, kernel, iterations=1)
    
    floodSeed = np.bitwise_not(removeDirt.astype(np.bool_))
    floodSeed[1:-1,1:-1]=False
    Edgeonly = floodBySeed(removeDirt,floodSeed,0,1,0,growthMin = 0, growthMax=100000)
    Edgeonly = np.bitwise_not(Edgeonly).astype(np.uint8)
    #Edgeonly = cv2.morphologyEx(threshed, cv2.MORPH_OPEN, kernel, iterations=1)
    kernel2 = np.ones((scaledMaxFeat/2,scaledMaxFeat/2))
    Edgeonly = cv2.morphologyEx(Edgeonly, cv2.MORPH_ERODE, kernel2)
    EdgeMask = misc.imresize(Edgeonly,image.shape,interp='bilinear')
    #EdgeMask=EdgeMask.astype(image.dtype)
    return EdgeMask

#Test this out on a few others'''