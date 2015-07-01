# This module contains all the functions necessary for analyzing the dirt on a foil

import os
import cv2
import numpy as np

import GenSIP.functions as fun
import GenSIP.gencsv as gencsv
import GenSIP.measure as meas

####################################################################################

####################################################################################

def analyzeDirt(sss, res):
    """
     This is the master script that takes in a Sample Set String and runs the 
     comparison functions on all of the foil pictures in the before and after 
     folders of that sample set. 
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
    
    Data = {}
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
                (numbf,
                 numaf,
                 areabf,
                 areaaf,
                 perDirtLoss) = dirtVals
                (threshedbf, 
                 threshedaf) = dirtPicts
                                                
                # Create the output images
                # Save the threshed images for analysis
                cv2.imwrite('Output/Output_'+sss+'/DirtMaps/'+foilnum+' before_threshed.png',\
                threshedbf, [cv2.cv.CV_IMWRITE_PNG_COMPRESSION,6])
                cv2.imwrite('Output/Output_'+sss+'/DirtMaps/'+foilnum+' after_threshed.png',\
                threshedaf, [cv2.cv.CV_IMWRITE_PNG_COMPRESSION,6])
                
        	# Write the output to the Data Dictionary
        	Data[foilnum] = {'Dirt Count Before':numbf, 
        	                 'Dirt Count After':numaf, 
                                 'Dirt Area Before':areabf, 
                                 'Dirt Area After':areaaf, 
                                 'Approx Percent Dirt Loss by Area':perDirtLoss}
                                 
            else: print "No after image for foil "+ foilnum
        else:
	   print "Not a picture: " + fn
	   
	   
    """Write the results to a .csv file"""
    # Make dirt output csv file
    filePath = 'Output/Output_'+sss+'/dirt_output_'+sss+'.csv'
    dirtCSV = gencsv.DataToCSV(filePath, sss)
    
    # Make column titles:
    ColHeaders = ['Foil #',
                  'Dirt Count Before', 
                  'Dirt Count After', 
                  'Dirt Area Before (mm^2)', 
                  'Dirt Area After (mm^2)', 
                  'Approx Percent Dirt Loss by Area']
                  
    # Write the Data to the CSV file
    dirtCSV.writeDataFromDict(Data,FirstColHead='Foil #',colHeads=ColHeaders)
    dirtCSV.closeCSVFile()
    
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
			
def dirtComp (beforeImg, afterImg, res, MaskEdges=True):
    """
    This is the dirt compare foils function. It takes the before and after images and 
    returns the number of dirt particles and area of dirt on each foil, and 
    the thresholded images used to calculate the amount of dirt.   
    """
    
    # Dirt analysis
    (numBef, 
     areaBef, 
     threshedBef, 
     sizesBef) = dirtnalysis (beforeImg, 
                              res, 
                              MaskEdges=True, 
                              retSizes=True)
    (numAft, 
     areaAft, 
     threshedAft, 
     sizesAft) = dirtnalysis (afterImg, 
                              res, 
                              MaskEdges=True, 
                              retSizes=True)
                              
    # Calculate difference in area
    if areaBef!=0: # Prevents division by zero
        perDirtLoss = round(100*(areaBef-areaAft)/float(areaBef),1)
    else:
        perDirtLoss = np.nan
    
    # put all results into return tuples
    ret = (numBef,
           numAft,
           areaBef,
           areaAft, 
           perDirtLoss)
           
    picts = (threshedBef, 
             threshedAft)
             
    return ret, picts

####################################################################################

####################################################################################

def dirtnalysis (img, res, MaskEdges=True, retSizes=False):
    """
    Performs dirt analysis on one image. 
    """
    
    # Dirt analysis
    threshed, masked = isolateDirt(img)
    area,num,sizes,labelled = meas.calcDirt(threshed, 
                                            res, 
                                            returnSizes=True,
                                            returnLabelled=True,
                                            getAreaInSquaremm=True)
    area = round(area,5)
    threshed = (masked/255)*np.bitwise_not(threshed)

    # put all results into return tuples
    if retSizes:
        return num, area, threshed, sizes
    else:
        return num, area, threshed

####################################################################################

####################################################################################

def isolateDirt (img, MaskEdges=True):
    poster = fun.makePoster(img)
    threshed,masked = fun.regionalThresh(img, poster,
                                         p=8, 
                                         d=28,
                                         m=55, 
                                         pt=60,
                                         MaskEdges=MaskEdges,
                                         returnMask=MaskEdges,
                                         MoDirt='dirt')
    return threshed, masked
'''
EXTRA and OUTDATED CODE:
    
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
	

'''