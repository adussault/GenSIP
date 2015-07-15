# -*- coding: utf-8 -*-
"""
This module contains the functions necessary for analyzing the dirt on a foil used
in the cleantests method. 

"""
import os
import cv2
import numpy as np

import GenSIP.functions as fun
import GenSIP.gencsv as gencsv
import GenSIP.measure as meas

####################################################################################

####################################################################################

def analyzeDirt(sss, res, verbose=False):
    """
     This is the master script that takes in a Sample Set String and runs the 
     comparison functions on all of the foil pictures in the before and after 
     folders of that sample set. 
        Inputs:
        - sss -  The Sample Set String, a short identifier for whichever set
                of SEM scans you are running.
        - res - Resolution of the image, in square microns per pixel. 
        Key-Word Arguments:
        - verbose = False - prints verbose output if set to True.
   
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
                if verbose:
                    print "Now comparing before and after for foil "+str(foilnum)
                # Compare the two images:
                dirtVals, dirtPicts, dirtSizes = dirtComp(beforeImg, afterImg, 
                                                          res, retSizeData=True,
                                                          verbose=verbose)   
                (numbf,
                 numaf,
                 areabf,
                 areaaf,
                 perDirtLoss) = dirtVals
                 
                (threshedbf, 
                 threshedaf) = dirtPicts
                 
                (BefMean, 
                 AftMean,
                 BefMax, 
                 AftMax,
                 BefOver100, 
                 AftOver100) = dirtSizes
                                                            
                # Create the output images
                # Save the threshed images for analysis
                cv2.imwrite('Output/Output_'+sss+'/DirtMaps/'+foilnum+' before_threshed.png',\
                threshedbf, [cv2.cv.CV_IMWRITE_PNG_COMPRESSION,6])
                cv2.imwrite('Output/Output_'+sss+'/DirtMaps/'+foilnum+' after_threshed.png',\
                threshedaf, [cv2.cv.CV_IMWRITE_PNG_COMPRESSION,6])
                
        	# Write the output to the Data Dictionary
        	Data[foilnum] = {'Dirt Count Before':numbf, 
        	                 'Dirt Count After':numaf, 
                                 'Dirt Area Before (mm^2)':areabf, 
                                 'Dirt Area After (mm^2)':areaaf, 
                                 'Approx % Dirt Loss by Area':perDirtLoss,
                                 'Mean Part. Area Before (micron^2)':BefMean,
                                 'Mean Part. Area After (micron^2)':AftMean,
                                 'Max Part. Area Before (micron^2)':BefMax,
                                 'Max Part. Area After (micron^2)':AftMax,
                                 'Approx % Parts. w/ >100micron diam. Before':BefOver100,
                                 'Approx % Parts. w/ >100micron diam. After':AftOver100
                                 }
                                 
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
                  'Approx % Dirt Loss by Area',
                  'Mean Part. Area Before (micron^2)',
                  'Mean Part. Area After (micron^2)',
                  'Max Part. Area Before (micron^2)',
                  'Max Part. Area After (micron^2)',
                  'Approx % Parts. w/ >100micron diam. Before',
                  'Approx % Parts. w/ >100micron diam. After']
                  
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
			
def dirtComp (beforeImg, afterImg, res, MaskEdges=True, retSizeData=False, verbose=False):
    """
    This is the dirt compare foils function. It takes the before and after images and 
    returns the number of dirt particles and area of dirt on each foil, and 
    the thresholded images used to calculate the amount of dirt. 
    
    Inputs: 
        - img - image as a numpy.ndarray
        - res - resolution of the image in square microns/pixel
    Key-Word Arguments:
        - MaskEdges = True - option to automatically mask off the background if 
            set to True
        - retSizeData = False - option to return the dirt size data if set to True
        - verbose = False - prints verbose output if set to True.
    Returns 2 tuples (or 3 if retSizeData = True):
        • ret - numerical dirt comparison data:
            - numBef - number of dirt particles
            - numAft
            - areaBef - total area of dirt particles in square microns
            - areaAft
        • picts - the thresholded images
            - threshedBef
            - threshedAft
        • sizeData - 
            - BefMean - mean dirt particle size in square microns
            - AftMean
            - BefMax - maximum dirt particle size in square microns
            - AftMax
            - BefPercOver100 - Percent of particles with an approximate diameter 
                               greater than 100 microns
            - AftPercOver100
    """
    
    # Dirt analysis
    (numBef, 
     areaBef, 
     threshedBef, 
     sizesBef) = dirtnalysis (beforeImg, 
                              res, 
                              MaskEdges=True, 
                              retSizes=True,
                              verbose=verbose)
    (numAft, 
     areaAft, 
     threshedAft, 
     sizesAft) = dirtnalysis (afterImg, 
                              res, 
                              MaskEdges=True, 
                              retSizes=True,
                              verbose=verbose)
                              
    BefMean, BefMax, BefPercOver100 = meas.getDirtSizeData(sizesBef, res)
    AftMean, AftMax, AftPercOver100 = meas.getDirtSizeData(sizesAft, res)
    
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
    
    sizeData = (BefMean, 
                AftMean,
                BefMax, 
                AftMax,
                BefPercOver100, 
                AftPercOver100)
                
    if retSizeData:
        return ret, picts, sizeData
    else:
        return ret, picts

####################################################################################

####################################################################################

def dirtnalysis (img, res, MaskEdges=True, retSizes=False, verbose=False):
    """
    Runs molybdenum analysis on a given image. 
        Inputs: 
        - img - image as a numpy.ndarray
        - res - resolution of the image in square microns/pixel
        Key-Word Arguments:
        - MaskEdges = True - option to automatically mask off the background if 
            set to True
        - retSizes = False - option to returnt the dirt size data if set to True
        - verbose = False - prints verbose output if set to True.
        Returns a tuple containing:
            num - number of dirt particles
            area - area of the dirt in the image in square microns
            threshed - the dirt thresholded image (white dirt on black background) 
                as a numpy ndarray
            sizes[optional] -  a 1-dimensional numpy ndarray listing out the sizes
               (area) of each dirt particle in pixels 
    """
    
    # Dirt analysis
    threshed, masked = isolateDirt(img, verbose=verbose)
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

def isolateDirt (img, MaskEdges=True, verbose=False):
    """Performs the thresholding on an image"""
    poster = fun.makePoster(img)
    threshed,masked = fun.regionalThresh(img, poster,
                                         p=8, 
                                         d=28,
                                         m=55, 
                                         pt=60,
                                         MaskEdges=MaskEdges,
                                         returnMask=MaskEdges,
                                         MoDirt='dirt',
                                         verbose=verbose)
    return threshed, masked
