"""
# GenTest is a module containing functions for testing out the functions in 
# GenSIP on a given folder of test images. 
Abandoned.
"""

import os
import cv2
import numpy as np
from matplotlib import pyplot as plt

direct = "FoilPicts/After/After_TEST/"
def original(dirpath):
# This is a complementary function to testRun. It simply runs a function on a  
#  folder of images and returns the result. It is designed to compare the original
#  operation on an image with an adjusted operation or an operation with addi-
#  -tional steps.
    images = {}
    names = []
    for imagepathname in os.listdir(dirpath):
            image = cv2.imread(dirpath+imagepathname, cv2.CV_LOAD_IMAGE_GRAYSCALE)
            name = os.path.splitext(imagepathname)[0]
            names.append(name) 
            #test a function on the image
            result = maskEdge(image)
            result = filterAndThresh(result)
            #store edited image and name in dictionary
            images[name] = result
    return images, names
    
def testRun(dirpath):
# This function tests out a given image processing function on a set of test 
#  images in a given directory. Stores names of images and results of applying 
#  the function to the images in a dictionary. Returns the dictionary along with 
#  a corresponding list of picture names for easier access to the dictionary. 
    images = {}
    names = []
    for imagepathname in os.listdir(dirpath):
        image = cv2.imread(dirpath+imagepathname, cv2.CV_LOAD_IMAGE_GRAYSCALE)
        name = os.path.splitext(imagepathname)[0]
        names.append(name) 
        #Do some preprocessing to image
        masked = maskEdge(image)

        #test a function on the image
        result = filterAndThresh(masked)
        #store edited image and name in dictionary
        images[name] = result
    return images, names


def filterAndThresh(img, gaussBlur=3):
#Separate function for applying the proper filtering and thresholding to the cropped foil image so the dirt can be counted
# Image must be imported as grayscale.
# gaussBlur is the variable associated with the size of the Gaussian Blur, default set to 3
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


def dictshow(images, names, **kwargs):
# This is a useful function for plotting and displaying any number of images quickly 
#     color is the color option ("gray", "jet" etc.) used in plt.imshow
#     Rows is the number of rows of subplots of images
# Any number of images can be inserted
        # Assign kwargs
        rows = kwargs.get('rows',3)
        color = kwargs.get('color',"gray")
	numImgs = len(images)
	# Columns of images is equal to the number of images divided by the specified number 
	# of rows + the modulus of numImgs and rows.
	columns = numImgs/rows + numImgs%rows
	pos = range(len(images))
	for i in pos:
		plt.subplot(rows,columns,i+1),plt.imshow(images[names[i]],color)
	plt.show()
def funcShow(f):
    def new_func (*args,**kwargs):
        rtn = f(*args,**kwargs)
        show(rtn)
        return rtn
    return new_func
        
# This is a set of functions for timing other functions.
import time

def timedRun(func, *args, **kwargs):
    t1=time.time()
    func(*args,**kwargs)
    t2 = time.time()-t1
    convTime(t2)
    
def convTime(t):
    mins = round(t/60,0)
    seconds = round(t%60,3)
    print str(mins)+"mins "+str(seconds)+"secs"
    
def timed(f):
# This is a decorator that creats a timed version of whatever function is passed
#  through it. 
    def new_func (*args, **kwargs):
        t1=time.time()
        rtn = f(*args, **kwargs)
        t2 = time.time()-t1
        convTime(t2)
        return rtn
    return new_func
    
#def clearOutput (sss):
    #Clears the output of a set of images
