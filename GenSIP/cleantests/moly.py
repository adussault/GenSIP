# This module contains all the functions necessary for analyzing the molybdenum loss
# on the foils. 

import cv2
import os
import numpy as np

import GenSIP.functions as fun
import GenSIP.measure as meas
import GenSIP.gencsv as gencsv

####################################################################################

####################################################################################

def analyzeMoly (sss, res, verbose=False):
    """ 
    This is the master function for comparing the exposed Pt on foil images.
    It takes in a 'Sample Set String' ('sss') and runs the comparison functions 
    on all of the foil pictures in the before and after folders of that sample set. 
    res is the resolution of the image, as micron per pixel. default set to one.
        Inputs:
        - sss -  The Sample Set String, a short identifier for whichever set
                of SEM scans you are running.
        - res - Resolution of the image, in square microns per pixel.     
    """
        
    # Declare a list of acceptable file types: tif and jpg
    filetypes = ['.tif', '.jpg', '.jpeg','.tiff']
    # Make list of files in the corresponding Before and After folders
    befpics = sorted(os.listdir('InputPicts/Before/Before_'+sss))
    aftpics = sorted(os.listdir('InputPicts/After/After_'+sss))
    
    # Make output folder if does not exist
    if not os.path.exists('Output/Output_'+sss):
        os.makedirs('Output/Output_'+sss)
    if not os.path.exists('Output/Output_'+sss+'/PtMaps'):
        os.makedirs('Output/Output_'+sss+'/PtMaps')
        
    Data = {}
    # fn stands for filename. befpics is a list of the filenames in the folder 
    # containing the before pictures
    for fn in befpics:
        filename,exten = os.path.splitext(fn)
	befpath = 'InputPicts/Before/Before_'+sss+'/'+fn
	if exten in filetypes:
	    # Get number of current foil
            foilnum = filename.split()[0]
            if verbose: print "Now comparing foil %s" %foilnum
            aftfoil = foilnum + ' after clean' + exten
            aftpath = 'InputPicts/After/After_'+sss+'/'+aftfoil
            # Check to see if the corresponding picture is in the after folder
            if aftfoil in aftpics:
                befImg = fun.loadImg(befpath, cv2.CV_LOAD_IMAGE_GRAYSCALE)
                aftImg = fun.loadImg(aftpath, cv2.CV_LOAD_IMAGE_GRAYSCALE)
                # Compare the two images:
                Ptvals, PtPicts = MoComp(befImg,aftImg, res,verbose=verbose)
		(PtAreabf,
		 PtAreaaf,
		 areaLoss,
		 Moloss,
		 PctMo) = Ptvals
		(Ptbef, 
		 Ptaft) = PtPicts
		
		# Create the output images
		cv2.imwrite('Output/Output_'+sss+'/PtMaps/'+foilnum+' before_threshed.png',\
		 Ptbef, [cv2.cv.CV_IMWRITE_PNG_COMPRESSION,6])
		cv2.imwrite('Output/Output_'+sss+'/PtMaps/'+foilnum+' after_threshed.png',\
		 Ptaft, [cv2.cv.CV_IMWRITE_PNG_COMPRESSION,6])
		 
		# Write the output to a new line in the csv file
		Data[foilnum] = {'Pt Area Before (mm^2)':PtAreabf, 
		                 'Pt Area After (mm^2)':PtAreaaf, 
                                 'Area of Mo Loss (mm^2)':areaLoss, 
                                 'Approx Mo Loss (micrograms)':round(Moloss,2),
                                 '% Mo lost':PctMo}
            else: print "No after image for foil "+foilnum
        else:
            print "Not a picture: %s" % filename
            
    """Write results Data to a .csv file"""
    # Make dirt and Pt output csv files and corresponding writers
    filePath = 'Output/Output_'+sss+'/Mo_output_'+sss+'.csv'
    PtCSV = gencsv.DataToCSV(filePath, sss)
    
    # Make column titles:
    ColHeaders = ['Foil #',
                  'Pt Area Before (mm^2)', 
                  'Pt Area After (mm^2)', 
                  'Area of Mo Loss (mm^2)', 
                  'Approx Mo Loss (micrograms)',
                  '% Mo lost']
                  
    # Write the Data to the CSV file
    PtCSV.writeDataFromDict(Data,colHeads=ColHeaders)
    PtCSV.closeCSVFile()
    
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
"""
def analyzeMoly_Single(sss,res,verbose=False):
    # Declare a list of acceptable file types: tif and jpg
    filetypes = ['.tif', '.jpg', '.jpeg','.tiff']
    # Make list of files in the corresponding Before and After folders
    befpics = sorted(os.listdir('InputPicts/Before/Before_'+sss))
    aftpics = sorted(os.listdir('InputPicts/After/After_'+sss))
    
    # Make output folder if does not exist
    if not os.path.exists('Output/Output_'+sss):
        os.makedirs('Output/Output_'+sss)
    if not os.path.exists('Output/Output_'+sss+'/PtMaps'):
        os.makedirs('Output/Output_'+sss+'/PtMaps')
        
    Data = {}
    # fn stands for filename. befpics is a list of the filenames in the folder 
    # containing the before pictures
    for fn in befpics:
        filename,exten = os.path.splitext(fn)
	befpath = 'InputPicts/Before/Before_'+sss+'/'+fn
	if exten in filetypes:
	    # Get number of current foil
            foilnum = filename.split()[0]
            if verbose: print "Now comparing foil %s" %foilnum
            aftfoil = foilnum + ' after clean' + exten
            aftpath = 'InputPicts/After/After_'+sss+'/'+aftfoil
            # Check to see if the corresponding picture is in the after folder
            if aftfoil in aftpics:
                befImg = fun.loadImg(befpath, cv2.CV_LOAD_IMAGE_GRAYSCALE)
                aftImg = fun.loadImg(aftpath, cv2.CV_LOAD_IMAGE_GRAYSCALE)
                # Compare the two images:
                Ptvals, PtPicts = MoComp(befImg,aftImg, res,verbose=verbose)
		(PtAreabf,
		 PtAreaaf,
		 areaLoss,
		 Moloss,
		 PctMo) = Ptvals
		(Ptbef, 
		 Ptaft) = PtPicts
		
		# Create the output images
		cv2.imwrite('Output/Output_'+sss+'/PtMaps/'+foilnum+' before_threshed.png',\
		 Ptbef, [cv2.cv.CV_IMWRITE_PNG_COMPRESSION,6])
		cv2.imwrite('Output/Output_'+sss+'/PtMaps/'+foilnum+' after_threshed.png',\
		 Ptaft, [cv2.cv.CV_IMWRITE_PNG_COMPRESSION,6])
		 
		# Write the output to a new line in the csv file
		Data[foilnum] = {'Pt Area Before (mm^2)':PtAreabf, 
		                 'Pt Area After (mm^2)':PtAreaaf, 
                                 'Area of Mo Loss (mm^2)':areaLoss, 
                                 'Approx Mo Loss (micrograms)':round(Moloss,2),
                                 '% Mo lost':PctMo}
            else: print "No after image for foil "+foilnum
        else:
            print "Not a picture: %s" % filename
            
    '''Write results Data to a .csv file'''
    # Make dirt and Pt output csv files and corresponding writers
    filePath = 'Output/Output_'+sss+'/Mo_output_'+sss+'.csv'
    PtCSV = gencsv.DataToCSV(filePath, sss)
    
    # Make column titles:
    ColHeaders = ['Foil #',
                  'Pt Area Before (mm^2)', 
                  'Pt Area After (mm^2)', 
                  'Area of Mo Loss (mm^2)', 
                  'Approx Mo Loss (micrograms)',
                  '% Mo lost']
                  
    # Write the Data to the CSV file
    PtCSV.writeDataFromDict(Data,colHeads=ColHeaders)
    PtCSV.closeCSVFile()
    
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
"""


####################################################################################

####################################################################################

def MoComp (beforeImg, afterImg, res, verbose=False):
    """
    This is the exposed platinum compare foils function. It takes the pathnamess of two images (before
    and after) and the image resolution, in square microns per pixel. Default set to 1.
    NOTE: This function assumes before and after images have the same resolution!
    MoComp returns two tuples: 
       -The first tuple contains the area of exposed platinum before and after 
        cleaning, the area of moly loss, the micrograms of moly lost, and the 
        approximate % of the original molybdenum lost.
       -The second tuple contains the binary before and after images of the
        exposed platinum used in these calculations. 
    """

    BEFDATA = Monalysis(beforeImg,res,verbose=verbose)
    (PtAreaBef, 
    FoilAreaBef, 
    MolyAreaBef,
    MolyMassBef, 
    PtImgBef) = BEFDATA
    
    AFTDATA = Monalysis(afterImg,res,verbose=verbose)
    (PtAreaAft, 
    FoilAreaAft, 
    MolyAreaAft,
    MolyMassAft, 
    PtImgAft) = AFTDATA
    
    # Calculate difference in area of exposed Pt:
    # Normalize values to the area of the before image
    befaftRatio = FoilAreaBef/FoilAreaAft
    areaLoss = PtAreaAft*befaftRatio-PtAreaBef
    if verbose:
        print "Total foil area before: %s" %str(FoilAreaBef)
        print "Total foil area after:  %s" %str(FoilAreaAft)
        print "Ratio of before to after: %s" %str(befaftRatio)
        print "Mo area before: %s" %str(MolyAreaBef)
        print "Mo area after:  %s" % str(MolyAreaAft)
    pcntMolyLoss = 1-(MolyAreaAft/float(MolyAreaBef))*befaftRatio
    pcntMolyLoss = round(pcntMolyLoss*100,1)
    # Now approximate molybdenum loss:
    # assuming 1pixel == 1 micron
    # Moly thickness ~300nm = .3micron, moly density = 10.2 g/cm^3 = 10.2e-9 mg/micron^3
    molyLoss = areaLoss*.3*10.2 #moly loss in micrograms
    molyLoss = round(molyLoss,2)
   	
    # Convert area values to mm^2 instead of microns^2, and round:
    PtAreaBef = round(PtAreaBef, 4)
    PtAreaAft = round(PtAreaAft, 4)
    areaLoss = round(areaLoss, 4)
    
    # put all results into return tuples
    ret = (PtAreaBef, 
           PtAreaAft, 
           areaLoss, 
           molyLoss, 
           pcntMolyLoss)
           
    picts = (PtImgBef, 
             PtImgAft)
             
    return ret, picts
    
####################################################################################

####################################################################################

def Monalysis(img, res, verbose=False):

    # Generate binary thresholds for Platinum:
    PtImg = isolatePt(img)

    # Approximate the percent of molybdenum lost:
    # Get the approximate foil area in square millimeters
    FoilArea = fun.getFoilArea(img, res, getAreaInSquaremm=True)

    # Calculate area of exposed platinum, adjust for resolution
    # Result is area in square millimeters
    PtArea = meas.calcExposedPt(PtImg, res)

    # Calculate difference in area of exposed Pt:
    # Normalize values to the area of the before image
    MolyArea = FoilArea-PtArea
    # Moly thickness ~300nm = .3micron, moly density = 10.2 g/cm^3 = 10.2e-9 mg/micron^3
    MolyMass = MolyArea*.3*10.2 #moly mass in micrograms
   	
    # Convert binary images to uint8 and make all ones into 255
    PtImg = PtImg.astype(np.uint8)
    PtImg = PtImg*255
    
    # put all results into return tuples
    ret = (PtArea, 
           FoilArea, 
           MolyArea, 
           MolyMass, 
           PtImg)
           
    return ret



####################################################################################

####################################################################################

def isolatePt (image):
    """
    This function filters and thresholds the image using regionalThres in order 
    to estimate the area of exposed Pt. Argument "image" must be ndarray.
    """
    poster = fun.makePoster(image)
    
    # Threshold the image. This is a global threshold. There is probably a better one out there.
    isoPt = fun.regionalThresh(image, poster,
                               p=90,
                               d=120,
                               m=180,
                               pt=180,
                               gaussBlur=3,
                               MaskEdges=False,
                               MoDirt='mo')
                               
    # Make the image into a boolean image
    isoPt = isoPt.astype(np.bool_)
    
    return isoPt

####################################################################################

####################################################################################
'''
EXTRA and OUTDATED CODE:
    	
def calcPtArea(img, comp=True):
    """
    This function calculates the area of Pt in a binary thresholded image of the foil.
    Calculate WHITE or BLACK area in a binary (boolean) image by counting pixels.
    """
    if str(img.dtype) == 'bool':
        if comp == False:
            inv = np.logical_not(img)
            return np.sum(inv)
        else:
            return np.sum(img)
    else:
        print "Only accept binary images!"

'''