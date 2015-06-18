# This module contains all the functions necessary for analyzing the molybdenum loss
# on the foils. 

import GenSIP.functions as fun
import GenSIP.measure as meas
import cv2
import os
import csv
import numpy as np
from socket import gethostname

####################################################################################

####################################################################################

def analyzeMoly (sss, res):
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
    
    # Make dirt and Pt output csv files and corresponding writers
    Ptdata = open('Output/Output_'+sss+'/Mo_output_'+sss+'.csv','w+b')
    Ptwriter = csv.writer(Ptdata)
	
    # Make header: 
    # First line is the title and the type of foil and test run:
    if sss.endswith("NF"):
        Ptwriter.writerow(["GenSIP Data","non-flight foils", '',sss.strip("NF")+" Test"])
    elif sss.endswith("F"):
        Ptwriter.writerow(["GenSIP Data","flight foils", '',sss.strip("F")+" Test"])
    else:
        Ptwriter.writerow(["GenSIP Data", "Molybdenum", sss])
    # Second row is the date and time of the run and the computer on which 
    # this function was run:
    datestring = fun.getDateString()
    host = os.path.splitext(gethostname())[0]
    Version = fun.getGenSIPVersion()
    Ptwriter.writerow(["Date:", datestring, "Computer:", host, "Version:", Version])
	
    # Make column titles:
    Ptwriter.writerow(['Foil #','Pt Area Before (mm^2)', 'Pt Area After (mm^2)', \
    'Area of Mo Loss (mm^2)', 'Approx Mo Loss (micrograms)','% Mo lost'])
    for fn in befpics:
        filename,exten = os.path.splitext(fn)
	befpath = 'InputPicts/Before/Before_'+sss+'/'+fn
	if exten in filetypes:
	    # Get number of current foil
            foilnum = filename.split()[0]
            print "Now comparing foil %s" %foilnum
            aftfoil = foilnum + ' after clean' + exten
            aftpath = 'InputPicts/After/After_'+sss+'/'+aftfoil
            # Check to see if the corresponding picture is in the after folder
            if aftfoil in aftpics:
                # Compare the two images:
                Ptvals, PtPicts =MoComp(befpath,aftpath, res)
		(PtAreabf,PtAreaaf,areaLoss,Moloss,PctMo) = Ptvals
		(Ptbef, Ptaft) = PtPicts

		# Convert binary images to uint8 and make all ones into 255
		Ptbef,Ptaft = Ptbef.astype(np.uint8),Ptaft.astype(np.uint8)
		Ptbef,Ptaft = Ptbef*255,Ptaft*255
		# Write the output to a new line in the csv file
		Ptwriter.writerow([foilnum, PtAreabf, PtAreaaf, areaLoss, round(Moloss,2),PctMo])
		# Create the output images
		# outfoil = compareMap(picts)
		# Save outfoil to Output folder
		# cv2.imwrite('Output/Output_'+sss+'/'+foilnum+' compare'+exten, outfoil)
		# Save the threshed images for analysis

		cv2.imwrite('Output/Output_'+sss+'/PtMaps/'+foilnum+' before_threshed.png',\
		 Ptbef, [cv2.cv.CV_IMWRITE_PNG_COMPRESSION,6])
		cv2.imwrite('Output/Output_'+sss+'/PtMaps/'+foilnum+' after_threshed.png',\
		 Ptaft, [cv2.cv.CV_IMWRITE_PNG_COMPRESSION,6])
            else: print "No after image for foil "+foilnum
        else:
            print "Not a picture: %s" % filename
			
####################################################################################

####################################################################################

def MoComp (before, after, res):
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

    # Load images
    bf = fun.loadImg(before, cv2.CV_LOAD_IMAGE_GRAYSCALE)
    af = fun.loadImg(after, cv2.CV_LOAD_IMAGE_GRAYSCALE)

    # Generate binary thresholds for Platinum:
    Ptbef,Ptaft = isolatePt(bf),isolatePt(af)

    # Approximate the percent of molybdenum lost:
    # Get the approximate foil area in square millimeters
    foilareabef = fun.getFoilArea(bf,res,getAreaInSquaremm=True)
    foilareaaft = fun.getFoilArea(af,res, getAreaInSquaremm=True)

    # Calculate area of exposed platinum, adjust for resolution
    # Result is area in square millimeters
    PtAreaBef = meas.calcExposedPt(Ptbef, res)
    PtAreaAft = meas.calcExposedPt(Ptaft, res)
	
    print "Total foil area before: %s" %str(foilareabef)
    print "Total foil area after:  %s" %str(foilareaaft)

    # Calculate difference in area of exposed Pt:
    # Normalize values to the area of the before image
    befaftRatio = foilareabef/foilareaaft
    areaLoss = PtAreaAft*befaftRatio-PtAreaBef
    MolyBef = foilareabef-PtAreaBef
    MolyAft = foilareaaft-PtAreaAft
    print "Mo area before: %s" %str(MolyBef)
    print "Mo area after:  %s" % str(MolyAft)
    pcntMolyLoss = 1-(MolyAft/MolyBef)*(foilareabef/foilareaaft)
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
    ret = (PtAreaBef, PtAreaAft, areaLoss, molyLoss,pcntMolyLoss)
    picts = (Ptbef, Ptaft)
    return ret, picts

####################################################################################

####################################################################################

def isolatePt (image):
    """
    This function filters and thresholds the image using regionalThres in order 
    to estimate the area of exposed Pt. Argument "image" must be ndarray.
    """
    poster = fun.makePoster(image)
    # Threshold the image. This is a global threshold. There is probably a better one out there.
    isoPt = fun.regionalThresh(image, poster,p=90,d=120,m=180,pt=180,gaussBlur=3,MoDirt='mo')
    # Make the image into a boolean image
    isoPt = isoPt.astype(np.bool_)
    return isoPt

####################################################################################

####################################################################################
	
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

