# This module contains all functions used across GenSIP
import cv2
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.figure as mplfig
from scipy import misc
from GenSIP.kuwahara import Kuwahara
import os
from time import localtime, asctime, struct_time

####################################################################################

####################################################################################

def maskEdge(img, thickness = 80):
    """
    This function masks off the outer edge of the foil, since the edge complicates dirt particle counting
    Parameters are:
        img = nparray of image,
        thickness = the number of pixels you want to take of the edge of the foil
                    Based on the image resolution, 1 micron ~= i pixel
    """
    #Load image
    image = img.copy()
    # Apply a morphological closing to the image to remove small dark areas:
    kernel = np.ones((3,3))
    cimage = cv2.morphologyEx(image, cv2.MORPH_CLOSE, kernel)
    #Threshold and blur image for contour detection
    blur2 = cv2.blur(cimage, (100,100))
    ret,thresh2 = cv2.threshold(blur2,50,255,cv2.THRESH_BINARY)
    # Put a black border around the thresholded image to prevent edge problems
    thresh2 = cv2.copyMakeBorder(thresh2,5,5,5,5,cv2.BORDER_CONSTANT,value=0)
    #Find the contours of the foil
    can1 = cv2.Canny(thresh2, 50, 200)
    contours,hierarchy = cv2.findContours(can1.copy(),cv2.RETR_TREE,cv2.CHAIN_APPROX_NONE)

    # Isolate and draw the largest contour, which should be the outline of the foil
    # Sort contours by size
    contours = sorted(contours, key = cv2.contourArea, reverse = True)
    # Largest contour should be the outer edge of foil
    outercnt = contours[0]
    #Create corresponding lists of contour areas and perimeters. 
    #cntAreas = list()
    #cntPerim = list()
    #Black out the edges of the foil
    #Double the mask_thickness to get the thickness of the masking line. 
    mask_thickness = thickness*2
    cv2.drawContours(image, [outercnt], -1, (0,255,0), mask_thickness)
    cv2.drawContours(thresh2, [outercnt], -1, (0,255,0), mask_thickness)
    # remove border from thresh2
    thresh2 = thresh2[5:thresh2.shape[0]-5,5:thresh2.shape[1]-5]
    # Return both the image and the foil area for the sake of future moly loss
    # approximations and for dirt counting.
    return image,thresh2

####################################################################################

####################################################################################
	
def makePoster(image,kern=6, KuSize=9,Gaus1=3,Gaus2=11,rsize=.1):
    """
    This method takes the image of the foil and creates a smoothed Kuwahara image
    used to make the poster for regional thresholding.
    """
    rsz = misc.imresize(image,rsize,interp='bicubic')
    gr = cv2.GaussianBlur(rsz, (Gaus1,Gaus1),0)
    kgr = Kuwahara(gr,KuSize)
    kgr = kgr.astype(np.uint8) # Make sure the Kuwahara image is uint8 so it doesn't scale
    rkgr = misc.imresize(kgr,(image.shape),interp='bicubic')
    grkgr = cv2.GaussianBlur(rkgr, (Gaus2,Gaus2),0)
    prkgr = posterfy(grkgr,kern)
    return prkgr

####################################################################################

####################################################################################

def posterfy(image,k_size=6):
    """
    Takes a gray image and sets all values within a given range to a single value
    this way we can divide up regions into primarily Pt, primarily Mo, dirt, or black
    """
    image_copy=image.copy()
    image_copy=image_copy.astype(np.uint8)
    blk=(image_copy<=4)
    pleat=(5<=image_copy)&(image_copy<=75)
    darkMo=(76<=image_copy)&(image_copy<=110)
    Mo=(111<=image_copy)&(image_copy<=199)
    Pt=(200<=image_copy)
    image_copy[blk]=0
    image_copy[pleat]=50
    image_copy[darkMo]=85
    image_copy[Mo]=150
    image_copy[Pt]=255
    #image_copy = maskEdge(image_copy,thickness=60)
    #Do a morphological opening step to eliminate the rough edges
    # I prefer opening over closing because it is more important to catch
    # all of the dark areas, as the Pt areas will have pretty sraightforward
    # threshold results
    kernel = np.ones((k_size,k_size))
    image_copy = cv2.morphologyEx(image_copy, cv2.MORPH_OPEN, kernel)
    return image_copy

####################################################################################

####################################################################################
    
def regionalThresh(ogimage,poster,p=8,d=28,m=55,pt=60,**kwargs):
    """
     This is the main thresholding method for analyzing the amount of Pt and dirt on
    the foils. It takes the posterized image and splits the image into regions based on
    the gray level in these different regions, then it assigns different threshold parameters 
    to the different regions, specified by the arguments p,d,m, and pt, which correspond to the
    regions 'pleat','dark Molybdenum', 'Molybdenum,' and 'platinum.'
    
     regionalThresh returns the final thresholded image with the black area around the foil
    cut out so that only the dirt appears. This allows mh.label to count the dirt and not get
    thrown off by the foil outline. If the option "GetMask" is set to True, then regionalThresh
    also returns the image of the outline of the foil and all regions that are black (<5).
    """
    gaussBlur=kwargs.get('gaussBlur',3)
    threshType=kwargs.get('threshType',0L)
    MaskEdges=kwargs.get('MaskEdges',0)
    returnMask=kwargs.get('returnMask',0)
    Mask=kwargs.get('Mask',0)
    MoDirt=kwargs.get('MoDirt','dirt')
    
    if poster.shape != ogimage.shape:
        raise Exception("The poster is not the same shape as the original image.")
        return
    Image = ogimage.astype(np.uint8)
    gPoster = poster.astype(np.uint8)
    threshedImage = np.zeros((Image.shape),dtype=np.uint8)

    # Create the images
    blk = gPoster.copy()
    pleat = gPoster.copy()
    darkMo = gPoster.copy()
    Mo = gPoster.copy()
    Pt = gPoster.copy()
    
    # isolate various sections
    blk[blk!=0]=1
    pleat[pleat!=50]=0
    darkMo[darkMo!=85]=0
    Mo[Mo!=150]=0
    Pt[Pt!=255]=0
    
    # mask the original image
    pleatMask = cv2.GaussianBlur(Image, (5,5), 0)
    pleatMask = (pleat/50)*pleatMask
    darkMoMask = cv2.GaussianBlur(Image, (5,5), 0)
    darkMoMask = (darkMo/85)*darkMoMask
    MoMask = cv2.GaussianBlur(Image, (gaussBlur,gaussBlur), 0)
    MoMask = (Mo/150)*MoMask
    PtMask = cv2.GaussianBlur(Image, (gaussBlur,gaussBlur), 0)
    PtMask = (Pt/255)*PtMask
    
    # apply adjusted threshold
    ret,pleat = cv2.threshold(pleatMask, p,255,threshType)
    ret,darkMo = cv2.threshold(darkMoMask, d,255,threshType)
    ret,MoMask = cv2.threshold(MoMask, m,255,threshType)
    ret,PtMask = cv2.threshold(PtMask, pt,255,threshType)
    # put it all together
    threshedImage = np.add(threshedImage, pleat)
    threshedImage = np.add(threshedImage, darkMo)
    threshedImage = np.add(threshedImage, MoMask)
    threshedImage = np.add(threshedImage, PtMask)
    
    # If using an external mask, apply it here:
    # Mask is assumed to be a binary image where black represents the areas to be
    # masked off. 
    if MaskEdges or type(Mask)==np.ndarray:
        
        if type(Mask)==np.ndarray:
            assert Mask.shape == ogimage.shape, \
            "Mask provided has different dimensions than the image to be threshed."
            threshMask = Mask.astype(np.bool_)*blk
            
        else:
            # If no mask is provided, use the maskEdge function
            masked,threshMask = maskEdge(ogimage)
            threshMask=threshMask.astype(np.bool_)*blk
        
        threshedImage = threshMask*threshedImage
        # Make sure the threshed image and mask are both np.uint8 images with 
        # white = 255.
        threshedImage[threshedImage!=0]=255
        threshMask = threshMask.astype(np.uint8)
        threshMask[threshMask!=0]=255
        
        if checkMoDirt(MoDirt)=='dirt':
            # If regionalThresh is working with dirt analysis, the inverse of the 
            # mask has to be applied so that all the masked off areas are whited
            # out and not blacked out. 
            threshedImage = threshedImage+np.bitwise_not(threshMask)
            # Then inverse the threshed image so it is white dirt on black:
            threshedImage = np.bitwise_not(threshedImage)      
        else:
            pass          
    else:
        if checkMoDirt(MoDirt)=='dirt':
            # Make sure the black areas do not show up as dirt. Inverse so white 
            # dirt is displayed over a black background.
            threshedImage = np.bitwise_not(threshedImage)*blk
            
    # Return image or tuple with image and mask if specified. 
    if returnMask:
        return threshedImage,threshMask
    else:
        return threshedImage
  
    ###############################
        # Combine the black region with the Mask
    """
    if type(Mask)==np.ndarray:
        threshMask = Mask.astype(np.bool_)
        masked = Image*threshMask
        comb = (blk.astype(np.bool_))*threshMask
        comb =comb.astype(np.uint8)
        comb[comb!=0]=255
        if MaskEdges:
            threshedImage = threshMask*threshedImage
            threshedImage = threshedImage+cv2.bitwise_not(comb)
    else:
        if MaskEdges:
        # Impose the masked edge:
            masked,threshMask = maskEdge(ogimage)
            threshMask[threshMask!=0]=1
            threshedImage = threshMask*threshedImage
        # Eliminate the black background so only dirt particles/Pt particles are visible
            comb = (blk/255)*threshMask
            threshedImage = threshedImage+cv2.bitwise_not(comb)
            
    if checkMoDirt(MoDirt) =='dirt':
        # Invert the image so it comes out as white dirt on black
        threshedImage = np.bitwise_not(threshedImage)
        # Dirt particle counting requires an 8-bit image with white as 255
        threshedImage = threshedImage.astype(np.uint8)
        threshedImage[threshedImage!=0]=255


    if returnMask:
        return threshedImage,comb
    else:
        return threshedImage
    """
####################################################################################

####################################################################################

def getFoilArea(img, res, getAreaInSquaremm=False):
    """
    FUNCTION MAY BE OUTDATED IF I REWORK MoComp AND isolatePT TO USE THE GetMask 
    OPTION AND USE THAT IMAGE TO APPROXIMATE FOIL AREA. 
    This is a function based on maskEdge that returns the approximate area of the foil.
        Inputs:
         - img = nparray of image,
         - res = the resolution of the image, in square microns per pixel
        Key-Word Arguments:
         - getAreaInSquaremm - If true, converts the getFoilArea from square microns
            to square millimeters.
    """
    #Load image
    image = img.copy()
    # Apply a morphological closing to the image to remove small dark areas:
    kernel = np.ones((3,3))
    cimage = cv2.morphologyEx(image, cv2.MORPH_CLOSE, kernel)
    #Threshold and blur image for contour detection
    blur2 = cv2.blur(cimage, (100,100))
    ret,thresh2 = cv2.threshold(blur2,50,255,cv2.THRESH_BINARY)
    thresh2 = cv2.copyMakeBorder(thresh2,5,5,5,5,cv2.BORDER_CONSTANT,value=0)
    # Approximate foil area is the sum of the white pixels divided by 255 and the 
    # resolution. This method of calculating the area has a weak point of being 
    # dependent on the white pixels, rather than the outer contour of the foil. 
    # This means that if parts of the foil are covered in enough dirt to make them
    # appear black, then they won't be added to the approximate foil area. 
    # For an area approximation that relies on the outer contour of the foil, 
    # see GenSIP.py in version 2 or previous, or in the "Outdated" folder of this 
    # version. That algorithm had the problem of occasionally messing up the outer
    # contour and throwing off the Moly loss calculations. 
    thresh2[thresh2!=0]=1
    foilarea = np.sum(thresh2)*res
    if getAreaInSquaremm:
        # Convert the foilarea to square mm
        foilarea = float(foilarea)*10**-6
    return foilarea

####################################################################################

####################################################################################

def show(*images, **kwargs):
    """
    This is a useful function for plotting and displaying any number of images quickly.
    Arguments:
        *images - any number of images you wish to display
    Key Word arguments:
        - color = 'gray' - designates the color option use in plt.imshow.
                Examples include "jet," "gray" (default), more in docstring
                for plt.imshow.
        - rows = 1 - designates numbers of rows to display the images on
    """
    # Get key word arguments:
    color = kwargs.get('color',"gray")
    rows = kwargs.get('rows',1)
    figure = kwargs.get('figure',False)
    numImgs = len(images)
    # Columns of images is equal to the number of images divided by the specified number 
    # of rows + the modulus of numImgs and rows.
    columns = numImgs/rows + numImgs%rows
    pos = range(len(images))
    if type(figure) == mplfig.Figure:
        subplts = {}
        for i in pos:
            subplts[i] = figure.add_subplot(rows,columns,i+1)
            subplts[i].imshow(images[i],color)
        figure.show()
    else:
        for i in pos:
       	    plt.subplot(rows,columns,i+1),plt.imshow(images[i],color)
        plt.show()

####################################################################################

####################################################################################


def loadImg (path, flag=cv2.CV_LOAD_IMAGE_GRAYSCALE):
    """
    This is a function for loading images. It basically just solves an issue with 
    cv2.imread and raises an exception if the path is wrong and cv2 returns a NoneType 
    rather than a numpy array.
    """
    image = cv2.imread(path, flag)
    if isinstance(image, type(None)):
    # check if image is an instance of type 'NoneType'
        if not(os.path.exists("path")):
        # Check if the path exists.
            raise Exception("Image file path does not exist.\n Path: {0}".format(path))
            return
        else:
        # This shouldn't come up often, but if the file exists but cv2.imread still
        # spits out a NoneType object, then this should catch the problem.
            raise Exception("Image file path exists, but cv2.imread couldn't load it.")
            return
    else:
        return image
        
####################################################################################

####################################################################################

def getDateString():
    """
    Get the date and time formated as "mm-dd-yyyy hh:minminAM"
    """
    t = localtime()
    if t.tm_hour < 12:
        end = "AM"
        hour = str(t.tm_hour)
    elif t.tm_hour>12:
        end = "PM"
        hour = str(t.tm_hour - 12)
    else:
        end = "PM"
        hour = str(t.tm_hour)
    datestring = str(t.tm_mon)+'-'+str(t.tm_mday)+'-'+str(t.tm_year)+' '+hour+'.'+str(t.tm_min)+end # Ex: "4-21-2015 3.50PM"
    return datestring
    
####################################################################################

####################################################################################

def justFrigginGetTheNormalTime(*args):
    """
    Returns time in the format: 'Mon Jun 22 17:49:39 2015'
    If called with no arguments, returns the current time in that format.
    If called with a time as a float or a struct_time object, returns that
    time in the specified format.
    """
    if args == []:
        retTime = asctime(localtime())
        return retTime
    else:
        retTimes = []
        for arg in args:
            if type(arg) == float:
                retTimes.append(asctime(localtime(arg)))
            elif type(arg) == struct_time:
                retTimes.append(asctime(arg))
        if len(retTimes) == 1:
            return retTimes[0]
        else:
            return retTimes
            
####################################################################################

####################################################################################

def getGenSIPVersion():
    """
    This returns the string that names the current version of GenSIP. It assumes 
    that the current working directory is the one containing the GenSIP folder
    """
    cwd = os.getcwd()
    directories = os.listdir(cwd)
    for f in directories:
        if f.startswith("GenSIP_v"):
            Version = f
            break
        elif f.startswith("GenSIP"):
            mod_time = justFrigginGetTheNormalTime((os.path.getmtime(f)))
            Version = mod_time + '('+f+')'
            break
    else:
        Version = "Unknown"
    return Version

####################################################################################

####################################################################################

def convTime(t):
    """
    Converts the time from seconds to a string in minutes and seconds
    """
    mins = round(t/60,0)
    seconds = round(t%60,3)
    return str(mins)+"mins "+str(seconds)+"secs"
    
####################################################################################

####################################################################################

def checkMoDirt(MoDirt):
    """
    Allows for several input options for MoDirt and converts them into a standard
    parameter 'mo' or 'dirt', and raises an exception if the input is invalid.
    """
    allowableMo = ("mo","moly","molybdenum","m")
    allowableDirt = ("dirt","d")
    if MoDirt.lower() in allowableMo:
        return 'mo'
    if MoDirt.lower() in allowableDirt:
        return 'dirt'
    else:
        raise Exception("\n MoDirt must be one of the following: \n"+\
        str(allowableMo)+"\n - - - - - or - - - - - \n"+str(allowableDirt))
    

####################################################################################

####################################################################################
def makeDiamondKernel(radius):
    kern = np.zeros((2*radius+1,2*radius+1))
    
    for i in range(radius+1):
        kern[radius-i:radius+i+1, i] = 1
        kern[i, radius+i+1:radius-i] = 1
        kern[-i, -radius-i:-radius-1+i] = 1
        kern[-radius-i:-radius-1+i, -i]=1
        #kern[:,radius-i]=1
    return kern.astype(np.uint8)
    
####################################################################################

####################################################################################
    
'''
##OUTDATED FUNCTIONS##
def filterAndThresh(img, gaussBlur=3):
#Separate function for applying the proper filtering and thresholding to the cropped foil image so the dirt can be counted
# Image must be imported as grayscale
# gaussBlur is the variabel associated with the size of the Gaussian Blur, default set to 3
	# Make a copy of the image so we do not affect the original. Not sure if this step is necessary.
	# imageF = filtered image
	imageF = img.copy()
	# Apply the filter to the image. A 3x3 Gaussian filter seems most appropriate as/
	# it allows for the detection of most small dirt particles while eliminating most/
	# anomalies associated with the Moly surface texture and scratches.
	# Will eliminate any dirt particles that are 2x1 pixels or less
	imageF = cv2.GaussianBlur(imageF, (gaussBlur,gaussBlur), 0)
	# Threshold the image
	# imageFT = filtered and thresholded image. 
	# Currently using Otsu's method for image thresholding
	ret,imageFT = cv2.threshold(imageF,0,255,cv2.THRESH_OTSU)
	# Return the filtered and thresholded image
	return imageFT
'''
def FizzBuzz():
    for i in range(1,101):
        if i%3==0:
            if i%5==0:
                print "FizzBuzz"
            else:
                print "Fizz"
        elif i%5==0:
            print "Buzz"
        else:
            print i