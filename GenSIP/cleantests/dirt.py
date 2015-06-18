# This module contains all the functions necessary for analyzing the dirt on a foil

import csv
import mahotas as mh
import os
import cv2
import numpy as np

import GenSIP.functions as fun
import GenSIP.measure as meas
from socket import gethostname

####################################################################################

####################################################################################

def analyzeDirt(sss, res):
    """
     This is the master script that takes in a Sample Set String and runs the comparison 
    functions on all of the foil pictures in the before and after folders of that sample set. 
        Inputs:
        - sss -  The Sample Set String, a short identifier for whichever set
                of SEM scans you are running.
        - res - Resolution of the image, in square microns per pixel. 
   
    """
    # Declare a list of acceptable file types: tif and jpg
    filetypes = ['.tif', '.jpg', '.jpeg','.tiff','.png','.bmp']
    # Make list of files in the corresponding Before and After folders
    befpics = sorted(os.listdir('InputPicts/Before/Before_'+sss))
    aftpics = sorted(os.listdir('InputPicts/After/After_'+sss))
    
    # Make output folder if does not exist
    if not os.path.exists('Output/Output_'+sss):
        os.makedirs('Output/Output_'+sss)
    if not os.path.exists('Output/Output_'+sss+'/DirtMaps'):
        os.makedirs('Output/Output_'+sss+'/DirtMaps')
    
    # Make dirt and Pt output csv files and corresponding writers
    dirtcsv = open('Output/Output_'+sss+'/dirt_output_'+sss+'.csv','w+b')
    dirtwriter = csv.writer(dirtcsv)
   	
        # Make header: 
    # First line is the title and the type of foil and test run:
    if sss.endswith("NF"):
        dirtwriter.writerow(["GenSIP Data","non-flight foils", '',sss.strip("NF")+" Test"])
    elif sss.endswith("F"):
        dirtwriter.writerow(["GenSIP Data","flight foils", '',sss.strip("F")+" Test"])
    else:
        dirtwriter.writerow(["GenSIP Data", "Dirt", sss])
	# Second row is the date and time of the run and the computer on which 
	# this function was run:
	datestring = fun.getDateString()
	host = os.path.splitext(gethostname())[0]
	Version = fun.getGenSIPVersion()
	dirtwriter.writerow(["Date:", datestring, "Computer:", host, "Version:", Version])
	
    # Make column titles:
    dirtwriter.writerow(['Foil #','Dirt Count Before', 'Dirt Count After', \
    'Dirt Area Before', 'Dirt Area After', 'Approx Percent Dirt Loss by Area'])

    # fn stands for filename. befpics is a list of the filenames in the folder 
    # containing the before pictures
    for fn in befpics:
        # Remove the extension from the filename, i.e. "40360,0202 before clean.tif"
        # becomes "40360,202 before clean" and ".tif"
        filename,exten = os.path.splitext(fn)
        befpath = 'InputPicts/Before/Before_'+sss+'/'+fn
        if exten in filetypes:
            # Get number of current foil by splitting the filename string and taking
            # the first element. I.e. "40360,0202 before clean" --> ["40360,0202","before", "clean"]
            # foilnum then is "40360,0202"
            foilnum = filename.split()[0]
            # Build the path to the after cleaning image
            aftfoil = foilnum + ' after clean' + exten
            aftpath = 'InputPicts/After/After_'+sss+'/'+aftfoil
            
            # Check to see if the corresponding picture is in the after folder
            if aftfoil in aftpics:
                # Load images
                beforeImg = cv2.imread(befpath, cv2.CV_LOAD_IMAGE_GRAYSCALE)
                afterImg = cv2.imread(aftpath, cv2.CV_LOAD_IMAGE_GRAYSCALE)
                # Compare the two images:
                dirtVals, dirtPicts = dirtComp(beforeImg,afterImg, res)   
                (numbf,numaf,areabf,areaaf) = dirtVals
                (threshedbf, threshedaf) = dirtPicts
                # Calculate difference in area
                if areabf!=0: # Prevents division by zero
                    perDirtLoss = 100*(areaaf-areabf)/areabf
        	else:
        	    perDirtLoss = np.nan
        	# Write the output to a new line in the csv file
        	dirtwriter.writerow([foilnum, numbf, numaf, areabf, areaaf, perDirtLoss])
                # Create the output images
                # Save the threshed images for analysis
                cv2.imwrite('Output/Output_'+sss+'/DirtMaps/'+foilnum+' before_threshed.png',\
                threshedbf, [cv2.cv.CV_IMWRITE_PNG_COMPRESSION,6])
                cv2.imwrite('Output/Output_'+sss+'/DirtMaps/'+foilnum+' after_threshed.png',\
                threshedaf, [cv2.cv.CV_IMWRITE_PNG_COMPRESSION,6])
            else: print "No after image for foil "+ foilnum
        else:
	   print "Not a picture: " + fn
	   
    # Check to see if any foils have an after picture but not a before picture:
    for fn in aftpics:
        filename,exten = os.path.splitext(fn)
        foilnum = filename.split()[0]
        beffoil = foilnum + ' before clean' + exten
        if exten in filetypes:
            if not(beffoil in(befpics)):
                print "No before image for foil "+ foilnum
        else:
            print "Not a picture: " + fn

####################################################################################

####################################################################################
			
def dirtComp (beforeImg, afterImg, res):
    """
    This is the dirt compare foils function. It takes the before and after images and 
    returns the number of dirt particles and area of dirt on each foil, and 
    the thresholded images used to calculate the amount of dirt.   
    """
    
    # Dirt analysis
    bposter = fun.makePoster(beforeImg)
    aposter = fun.makePoster(afterImg)
    bthreshed,bmasked = fun.regionalThresh(beforeImg,bposter,p=8,d=28,m=55,pt=60,MaskEdges=True,returnMask=True,MoDirt='dirt')
    athreshed,amasked = fun.regionalThresh(afterImg,aposter,p=8,d=28,m=55,pt=60,MaskEdges=True,returnMask=True,MoDirt='dirt')
    areabf,numbf,sizesbf,labbf = meas.calcDirt(bthreshed, res, returnSizes=True,returnLabelled=True)
    areaaf,numaf,sizesaf,labaf = meas.calcDirt(athreshed, res, returnSizes=True,returnLabelled=True)

    bthreshed = (bmasked/255)*np.bitwise_not(bthreshed)
    athreshed = (amasked/255)*np.bitwise_not(athreshed)

    # put all results into return tuples
    ret = (numbf,numaf,areabf,areaaf)
    picts = (bthreshed, athreshed)
    return ret, picts

####################################################################################

####################################################################################

def amountDirt(img, res):
    """
    Calculates the number of dirt particles and the area of the foil covered by dirt
    Takes a black and white image (black dirt on white foil)
    Returns the number of dirt particles, the area of dirt particles, and a labelled image
    inverse the colors since mahotas.label only works on white on black
    """
    inv = cv2.bitwise_not(img)
    #invDirt = cv2.bitwise_not(isoDirt(img,profile))
    labeledFoil,numDirt = mh.label(inv)
    # Don't count the foil in numDirt:
    numDirt -= 1
    #Calculate the area of the dirt using Findcontours
    sizes = mh.labeled.labeled_size(labeledFoil)
    # Sort sizes of particles by size in descending order:
    sizes = np.sort(sizes)[::-1]
    # Eliminate the foil from the "sizes" array:
    sizes = sizes[2:]
    # Total area of dirt is equal to the sum of the sizes.
    areaDirt = sum(sizes)*res
    return (numDirt, areaDirt, labeledFoil, sizes)
	

