"""
Contains functions that measure the exposed platinum, dirt particles, molybdenum,
etc. based on the thresholded images produced by analysis methods.
"""
import numpy as np
import cv2
import mahotas as mh
import os

import GenSIP.functions as fun
import matplotlib.pyplot as plt

####################################################################################

####################################################################################
             

def calcExposedPt (Ptimage, res,**kwargs):
    """
    Returns the area of exposed platinum in square millimeters.
    Input is the thresholded platinum image and the image resolution, 
    in square microns per pixel. 
    """
    getAreaInSquaremm=kwargs.get('getAreaInSquaremm',True)

    # Make sure the image minimum is 0 and the image is binary (all pixel values
    # are equal to either the maximum value or 0)
    if Ptimage.min()!=0 or Ptimage[(Ptimage!=0)&(Ptimage!=Ptimage.max())].size!=0:
    # Make sure the image is binary
        if Ptimage.max() == Ptimage.min():
            print "Warning: Image is all white! Platinum area set to zero."
            areaPt = 0
            return areaPt
        else:
            raise Exception(
            """
            ERROR: Image must be a binary image of 0 and a non-zero number.\n
            Image Max: {0} \n
            Image Min: {1} \n
            Other values: {2}""".format(Ptimage.max(),
                                        Ptimage.min(),
                                        Ptimage[(Ptimage!=Ptimage.min())&(Ptimage!=Ptimage.max())]))
    # Convert to uint8 format
    plat = Ptimage.astype(np.uint8)
    # Make sure the image is just made of 1s and 0s
    plat[plat!=0]=1
    # Area of Pt = sum of nonzero pixels x resolution x 10^-16
    areaPt = float(plat.sum())*res
    if getAreaInSquaremm: 
        areaPt = areaPt*10**-6
    return areaPt

####################################################################################

####################################################################################
    
def calcDirt(img, res, **kwargs):
    """
    Calculates the number of dirt particles and the area of the foil covered by dirt
    Takes a black and white image (white dirt on black foil)
    Returns the number of dirt particles, the area of dirt particles, and a labelled image
    inverse the colors since mahotas.label only works on white on black
    
        Key-word Arguments:
            returnSizes=False
                - Option that returns the sizes array as the third ouput if true
            returnLabelled = False
                - Option that returns the labelled array as the fourth output if true
            foilIncluded = False
                - If true than the input image includes the outline of the foil in
                black and the bakcground in white.
            areaInSqmm = False
                - Option to return the areaDirt in square mm rather than sq microns.
            BoundConds = np.ones((3,3))
                - Structuring element that tells the mahotas label function which
                nearest neighbors to consider as part of the same region. Default
                set to a 3x3 matrix of ones so it will consider the 8 nearest neighbors.
            minPartArea = 0
                - Minimum dirt particle area to be considered in the count and area
                approximation. In square microns. Default set to 0.
    
    """
    returnSizes=kwargs.get('returnSizes',False)
    returnLabelled=kwargs.get('returnLabelled',False)
    foilIncluded = kwargs.get('foilIncluded',False)
    getAreaInSquaremm=kwargs.get('getAreaInSquaremm',False)
    BoundConds = kwargs.get('BoundConds',np.ones((3,3)))
    minPartArea = kwargs.get('minPartArea',0)
    
    # Make sure the image is binary
    if img.max() == img.min():
        print "Warning: Image is all white! All dirt values set to zero."
        areaDirt=0
        numDirt=0
        sizes = np.zeros((0,)).astype(np.uint32)
        labeledFoil = np.zeros(img.shape)
            # Generate return tuple:
        ret = [areaDirt, numDirt]
        if returnSizes:
            ret.append(sizes)
        if returnLabelled:
            ret.append(labeledFoil)
        return tuple(ret)
        
    if img.min()!=0 or np.any(img[img!=img.min()]!=img.max()):
        raise Exception(
            """
            ERROR: Image must be a binary image of 0 and a non-zero number.\n
            Image Max: {0} \n
            Image Min: {1} \n
            Other values: {2}""".format(img.max(),
                                        img.min(),
                                        img[(img!=img.min())&(img!=img.max())]))
    
    #inv = cv2.bitwise_not(img)
    #invDirt = cv2.bitwise_not(isoDirt(img,profile))
    # Make a 3x3 matrix of ones as the structuring element so that any of the 8 nearest
    # neighbors are all considered part of the same region. 
    labeledFoil,numDirt = mh.label(img,Bc=BoundConds)
    
    #Calculate the area of the dirt using Findcontours
    sizes = mh.labeled.labeled_size(labeledFoil)
    # Sort sizes of particles by size in descending order:
    sizes = np.sort(sizes)[::-1]
    # Eliminate the background from the "sizes" array:
    # If the image includes the background and the foil, 
    # they will be largest and second-largest regions in the sizes area
    if foilIncluded:
        sizes = sizes[1:]
        # Don't count the foil in numDirt:
        numDirt -= 1
    # Otherwise just eliminate the background region:
    else:
        sizes = sizes[1:]
    # Consider the minimum in pixels:
    minPartPx = minPartArea/res
    sizes[sizes<minPartPx]=0
    sizes = sizes[sizes!=0]
    # Total area of dirt is equal to the sum of the sizes. 
    # In square microns unless otherwise specified.
    areaDirt = sum(sizes)*res
    # If specified, convert area to square mm
    if getAreaInSquaremm:
        areaDirt = float(areaDirt*10**-6)
    # Generate return tuple:
    ret = [areaDirt, numDirt]
    if returnSizes:
        ret.append(sizes)
    if returnLabelled:
        ret.append(labeledFoil)
    return tuple(ret)
    
####################################################################################

####################################################################################

def getDirtSizeData(DirtSizes, res):
    """
    Receives the dirt particle size array and the resolution of the image and 
    returns some data on the particle sizes.
    Inputs:
        - DirtSizes - 1-dimensional ndarray of dirt sizes in pixels.
        - res - resolution of the image in square microns/pixel area
    Returns a tuple containing:
        - MeanSize - Mean dirt particle area in square microns
        - MaxSize - Max dirt particle area in square microns
        - percAreaOver100 - Percent of particles with an approximate diameter 
                            greater than 100 microns
    """
    if DirtSizes.size==0:
        MeanSize = 0
        MaxSize = 0
        percAreaOver100 = 0
    else:
        DirtSizes = DirtSizes*res
        MeanSize = round(DirtSizes.mean(),1)
        MaxSize = DirtSizes.max()
        # Number of Particles with a diameter approximately over 100 microns,
        # Corresponding to an area of ~7854 square microns. 
        areaOver100 = DirtSizes[DirtSizes>7850].sum()
        percAreaOver100 = round(100*(float(areaOver100)/(DirtSizes.sum())),2)
        
    return MeanSize, MaxSize, percAreaOver100

####################################################################################

####################################################################################

def makeSizeHistogram(sizes, res, name, path):
    """
    Makes a histogram of the dirt size data
    Inputs:
        - sizes - a 1-dimensional numpy ndarray listing out the sizes (areas)
                    of each dirt particle in pixels 
        - res - resolution of the image in square microns/pixel area
        - name - name of the image to be the title of the histogram plot
        - path - path to folder to which the histogram plot will be saved. 
    """
    AreaSizes = sizes*res
    histoRange = AreaSizes.max()
    sizeHist,bins = np.histogram(AreaSizes,bins=histoRange)
    x = np.arange(1,histoRange+1,1)
    normHistoByArea = sizeHist*x
    FIG=plt.figure()
    PLT = FIG.add_subplot(111)
    FIG.suptitle(name+" Histograms")
    
    plt.xlabel("Particle Size in square Microns")
    plt.ylabel("Count")
    PLT.plot(x, normHistoByArea)
    FIG.savefig(os.path.join(path,name+"_DirtPartSize_Hist.png"))
    print "Histogram plot for "+name+" finished."
    plt.close()

####################################################################################

####################################################################################
               
# ERROR CLASSES:
class MeasError(Exception):
    """Generic error for the measure module"""
    pass
    
class AllWhiteError(MeasError):
    """Error if image is all white"""
    def __init__(self, msg):
        self.msg = msg
        
class AllBlackError(MeasError):
    """Error if image is all black"""
    def __init__(self, msg):
        self.msg = msg

    
              

                
                
                