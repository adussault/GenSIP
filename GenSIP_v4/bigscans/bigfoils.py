'''
 This is the main module for analyzing the SEM scans of the large foils. 
'''
import cv2
import os
import csv

import numpy as np
import mahotas as mh
import GenSIP_v4.functions as fun
import GenSIP_v4.measure as meas
import GenSIP_v4.bigscans.images as images


from socket import gethostname


#Q1 = fun.loadImg("InputPicts/FoilScans/Q1/panorama.tif",0)

################################################################################

################################################################################

def analyzePano(panPath, maskPath, res, foilname, Quarter="", MoDirt="Mo", GenPoster=False, verbose=True):
    """
    This function runs full analysis on a single panorama SEM scan of a foil. 
    Essentially all this does is split up the panorama and mask images, puts 
    the resulting subimages into folders, and then calls 'analyzeSubImages' on 
    the folders of sub-images. 
    
    The inputs are:
        Arguments:
            
        - panPath - the path string to the panorama image of the foil
        - maskPath - the path string to the panorama of a binary image denoting 
            areas to keep masked off from analysis. Currently the working proce
            -dure is for someone to manually mask off the edges and all cracks of
            the SEM scan in photoshop, and then save this mask panorama as a .png 
            or .tiff file. Hopefully in the future a more automated crack-and-
            -edge-finding algorithm can be developed.
        - res - Resolution of the image, in square microns per pixel. 
        - foilname - The name/serial number of the foil being analyzed, for 
            example: "40360,2". 
            
        Key-word Arguments:
        
        - Quarter - The section of the foil being analyzed, i.e. "Q1". 
            Default is set to an empty string. 
        - MoDirt - string designating whether to run dirt analysis or molybdenum
            analysis. Allowable inputs:
                "Mo","Moly","moly","molybdenum","Molybdenum","M","m"
                or for dirt analysis: "Dirt","dirt","D","d"
        - GenPoster - option if the user wants to generate and save poster images
            in order to see how the image is split up into regions. 
    """
    print "MoDirt:  " + MoDirt
    panorama = fun.loadImg(panPath, 0)
    images.splitImage(panorama, 256, path="InputPicts/FoilScans/"+foilname, extra=Quarter)
    del(panorama)
    
    mask = fun.loadImg(maskPath, 0)
    images.splitImage(mask, 256, path="InputPicts/FoilScans/"+foilname, extra=Quarter+"_mask")
    del(mask)
    
    panFolder = "InputPicts/FoilScans/"+foilname+"/sub_imgs_"+Quarter
    maskFolder = "InputPicts/FoilScans/"+foilname+"/sub_imgs_"+Quarter+"_mask"
    
    # Call analyze sub images. 
    analyzeSubImages(panFolder,maskFolder,res,foilname,Quarter,MoDirt,GenPoster)

################################################################################

################################################################################


def analyzeSubImages(panFolder,maskFolder,res, foilname,  Quarter="", MoDirt="Mo",  GenPoster=False, verbose= True):
    """
    This function runs the analysis on all of the subimages of the panorama. It 
    writes the output to a csv file and produces images of the dirt and exposed 
    platinum on the foil. 
    
        Arguments:
        - panFolder - Path string to the folder containing all of the sub-
            images of the panorama.
        - maskFolder - Path string to the folder containing all of the sub-
            images of the mask.
        - res - Resolution of the image, in square microns per pixel.
        - foilname - The name/serial number of the foil being analyzed, for 
            example: "40360,2". 
            
        Key-word Arguments:
        
        - Quarter - The section of the foil being analyzed, i.e. "Q1". 
            Default is set to an empty string. 
        - MoDirt - string designating whether to run dirt analysis or molybdenum
            analysis. Allowable inputs:
                "Mo","Moly","moly","molybdenum","Molybdenum","M","m"
                or for dirt analysis: "Dirt","dirt","D","d"
        - GenPoster - option if the user wants to generate and save poster images
            in order to see how the image is split up into regions. 
                
    """
    # Create a list of the the contents of the panFolder and maskFolder, which will 
    # Produce a list of the names of all of the subimages in the pan & mask folders.
    # I.e. one element of the list would be the string: "sub_004_004.tiff"
    
    panSubs = os.listdir(panFolder)
    maskSubs = os.listdir(maskFolder)
    
    # Make sure only subImages appear in panSubs and maskSubs
    panSubs = FILonlySubimages(panSubs, limitToType=0)
    maskSubs = FILonlySubimages(maskSubs, limitToType=0)
    
    # Create output folder if it does not exist
    outFolder = 'Output/Output_'+foilname+"_"+Quarter
    if not os.path.exists(outFolder):
        os.makedirs(outFolder)
    
    # Create output folder for the poster if it does not exist and is requested
    if GenPoster and not os.path.exists(outFolder+"/PosterMaps"):
        os.makedirs(outFolder+"/PosterMaps")
        
    ### MOLYBDENUM ________________________________________________
    ############################################################################
    if fun.checkMoDirt(MoDirt)=='mo':
        
        #Create output folder:
        if not os.path.exists(outFolder+'/PtMaps'):
            os.makedirs(outFolder+'/PtMaps')
        
        # Set up the output csv file:
        data = open(outFolder+'/'+Quarter+'_PtData.csv','w+b')
        Ptwriter = csv.writer(data)
        # Make Header
        Ptwriter.writerow(["GenSIP Pt Data","SEM Scan", foilname, Quarter])
        datestring = fun.getDateString()
        host = os.path.splitext(gethostname())[0]
        Version = fun.getGenSIPVersion()
        Ptwriter.writerow(["Date:", datestring])
        Ptwriter.writerow(["Computer:", host, "Version:", Version])
        Ptwriter.writerow(['SubImage #','Pt Area (mm^2)', 'Foil area','% Exposed Pt'])
        
        totFoilArea = 0
        totPtArea = 0
        
        for sub in panSubs:
            root,ext = os.path.splitext(sub)
            
            # Create the threshholded image
            subImage = fun.loadImg(panFolder+'/'+sub,0)
            subMask = fun.loadImg(maskFolder+'/'+sub,0)
            proc = images.bigPostPreProc(subImage)
            poster = images.bigPosterfy(proc)
            if GenPoster:
                cv2.imwrite(outFolder+'/PosterMaps/'+root+".png",poster, \
                [cv2.cv.CV_IMWRITE_PNG_COMPRESSION,6])
            threshed = bigRegionalThresh(subImage,poster,Mask=subMask,\
            p=150,d=180,m=210,hE=240,pt=253,MoDirt=MoDirt)
            threshed = threshed.astype(np.bool_)
            
            # Get amount of exposed Platinum and write output to csv file:
            
            PixPt = np.sum(threshed)
            AreaPt = round(PixPt*res*10**-6, 4)
            PixFoil = np.sum(subMask.astype(np.bool_))
            AreaFoil = round(PixFoil*res*10**-6, 4)
            
            if AreaFoil == 0:
                PercPt = np.nan
            else:
                PercPt = round(float(AreaPt)/float(AreaFoil)*100,2)
            Ptwriter.writerow([root,AreaPt,AreaFoil,PercPt])
            
            # Make output image
            threshed=threshed.astype(np.uint8)*255
            cv2.imwrite(outFolder+'/PtMaps/'+root+".png",threshed, \
            [cv2.cv.CV_IMWRITE_PNG_COMPRESSION,6])
            
            totPtArea = totPtArea + AreaPt
            totFoilArea = totFoilArea + AreaFoil
            if verbose: print sub
            del(threshed)
            
        totPtArea = round(float(totPtArea)/100,2)
        totFoilArea = round(float(totFoilArea)/100,2)
        if totFoilArea == 0:
            PercPt = np.nan
        else:
            PercPt = round(100*totPtArea/totFoilArea,2)
        Ptwriter.writerow(["TOTALS:"])
        Ptwriter.writerow(["Foil Area (cm^2)", totFoilArea])
        Ptwriter.writerow(["Exposed Pt Area (cm^2)", totPtArea])
        Ptwriter.writerow(["% Exposed Pt", PercPt])
        images.stitchImage(outFolder+'/PtMaps')
        
    ### Dirt ________________________________________________________________
    ############################################################################
    
    elif fun.checkMoDirt(MoDirt)=='dirt':
        
        #Create output folder:
        if not os.path.exists(outFolder+'/DirtMaps'):
            os.makedirs(outFolder+'/DirtMaps')
            
        # Set up the output csv file:
        data = open(outFolder+'/'+Quarter+'_dirtData.csv','w+b')
        dirtwriter = csv.writer(data)
        # Make Header
        dirtwriter.writerow(["GenSIP Dirt Data","SEM Scan", foilname, Quarter])
        datestring = fun.getDateString()
        host = os.path.splitext(gethostname())[0]
        Version = fun.getGenSIPVersion()
        dirtwriter.writerow(["Date:", datestring])
        dirtwriter.writerow(["Computer:", host, "Version:", Version])
        dirtwriter.writerow(['SubImage #',"Dirt Count","Dirt Area (microns^2)", "Foil area (mm^2)","% Covered in dirt"])
        
        totFoilArea = 0
        totDirtArea = 0
        
        for sub in panSubs:
            root,ext = os.path.splitext(sub)
            
            # Create the threshholded image
            subImage = fun.loadImg(panFolder+'/'+sub,0)
            subMask = fun.loadImg(maskFolder+'/'+sub,0)
            proc = images.bigPostPreProc(subImage)
            poster = images.bigPosterfy(proc)
            threshed = bigRegionalThresh(subImage,poster,Mask=subMask,\
            p=8,d=28,m=55,hE=60,pt=70,MoDirt=MoDirt)
            threshed = threshed.astype(np.uint8)
            
            # Make output image
            cv2.imwrite(outFolder+'/DirtMaps/'+root+".png",threshed, \
            [cv2.cv.CV_IMWRITE_PNG_COMPRESSION,6])
                        
            if GenPoster:
                cv2.imwrite(outFolder+'/PosterMaps/'+root+".png",poster, \
                [cv2.cv.CV_IMWRITE_PNG_COMPRESSION,6])
            
            # Get amount of dirt and write output to csv file:
            AreaDirt, numDirt,sizes,labeled = meas.calcDirt(threshed,res, \
            returnSizes=True,returnLabelled=True, getAreaInSquaremm=False)
            
            AreaDirt = round(float(AreaDirt), 6)
            PixFoil = np.sum(subMask.astype(np.bool_))
            AreaFoil = round(float(PixFoil)*res, 8)
            if AreaFoil == 0:
                PercDirt = np.nan
            else:
                PercDirt = round(float(AreaDirt)/float(AreaFoil)*100,2)
                
            dirtwriter.writerow([root,numDirt, AreaDirt,round(float(AreaFoil)*10**-6,2),PercDirt])
            
            totDirtArea = totDirtArea + AreaDirt*10**-6
            totFoilArea = totFoilArea + AreaFoil*10**-6
            if verbose: print sub +" dirt count: " + str(numDirt)

            del(threshed)
            
        # Calculate Totals and write the last few rows of the csv file. 
        totDirtArea = round(float(totDirtArea)/100,4)
        totFoilArea = round(float(totFoilArea)/100,2)
        if totFoilArea == 0:
                PercDirt = np.nan
        else:
            PercDirt = round(100*totDirtArea/totFoilArea,2)
        dirtwriter.writerow(["TOTALS:"])
        dirtwriter.writerow(["Foil Area (cm^2)", totFoilArea])
        dirtwriter.writerow(["Exposed Dirt Area (cm^2)", totDirtArea])
        dirtwriter.writerow(["% Covered in Dirt", PercDirt])
        data.close()
        images.stitchImage(outFolder+'/DirtMaps')
    if GenPoster:
        images.stitchImage(outFolder+'/PosterMaps')

################################################################################

################################################################################


def bigRegionalThresh(ogimage,poster,p=8,d=28,m=55,hE=60,pt=70,gaussBlur=3,threshType=0L,Mask=0,GetMask=0,MoDirt="Mo"):
    """
     This is the main thresholding method for analyzing the amount of Pt and dirt on
    the foils. It takes the posterized image and splits the image into regions based on
    the gray level in these different regions, then it assigns different threshold parameters 
    to the different regions, specified by the arguments p,d,m, and pt, which correspond to the
    regions 'pleat','dark Molybdenum', 'Molybdenum,' and 'platinum.' For the large scans,
    an additional region called "highEx" or high exposure, was created above Pt. 
    
     regionalThresh returns the final thresholded image with the black area around the foil
    cut out so that only the dirt appears. This allows mh.label to count the dirt and not get
    thrown off by the foil outline. If the option "GetMask" is set to True, then regionalThresh
    also returns the image of the outline of the foil and all regions that are black (<5).
    """
    if poster.shape != ogimage.shape:
        raise Exception("The two arrays are not the same shape.")
        return
    Image = ogimage.astype(np.uint8)
    gPoster = poster.astype(np.uint8)
    threshedImage = np.zeros((Image.shape),dtype=np.uint8)

    # Create the images
    blk = gPoster.copy()
    pleat = gPoster.copy()
    darkMo = gPoster.copy()
    Mo = gPoster.copy()
    highEx = gPoster.copy()
    Pt = gPoster.copy()

    # isolate various sections
    blk[blk!=0]=255
    pleat[pleat!=50]=0
    darkMo[darkMo!=85]=0
    Mo[Mo!=150]=0
    highEx[highEx!=200]=0
    Pt[Pt!=255]=0

    # mask the original image
    pleatMask = cv2.GaussianBlur(Image, (5,5), 0)
    pleatMask = (pleat/50)*pleatMask
    darkMoMask = cv2.GaussianBlur(Image, (5,5), 0)
    darkMoMask = (darkMo/85)*darkMoMask
    MoMask = cv2.GaussianBlur(Image, (gaussBlur,gaussBlur), 0)
    MoMask = (Mo/150)*MoMask
    highExMask = cv2.GaussianBlur(Image, (gaussBlur,gaussBlur), 0)
    highExMask = (highEx/200)*highExMask
    PtMask = cv2.GaussianBlur(Image, (gaussBlur,gaussBlur), 0)
    PtMask = (Pt/255)*PtMask

    # apply adjusted threshold
    ret,pleatMask = cv2.threshold(pleatMask, p,255,threshType)
    ret,darkMoMask = cv2.threshold(darkMoMask, d,255,threshType)
    ret,MoMask = cv2.threshold(MoMask, m,255,threshType)
    ret,highExMask = cv2.threshold(highExMask, hE,255,threshType)
    ret,PtMask = cv2.threshold(PtMask, pt,255,threshType)

    # put it all together
    threshedImage = np.add(threshedImage, pleatMask)
    threshedImage = np.add(threshedImage, darkMoMask)
    threshedImage = np.add(threshedImage, MoMask)
    threshedImage = np.add(threshedImage, highExMask)
    threshedImage = np.add(threshedImage, PtMask)


    # Apply Mask if provided
    if type(Mask)==np.ndarray and Mask.shape==ogimage.shape:
        blur = cv2.GaussianBlur(Mask, (gaussBlur,gaussBlur), 0).astype(np.bool_)
        #invMask = np.bitwise_not(blur).astype(np.bool)
        if fun.checkMoDirt(MoDirt) =='dirt':
            # Invert the image so it comes out as white dirt on black
            threshedImage = np.bitwise_not(threshedImage)
            threshedImage = threshedImage*blur
            # Dirt particle counting requires an 8-bit image with white as 255
            threshedImage = threshedImage.astype(np.uint8)
            threshedImage[threshedImage!=0]=255
        elif fun.checkMoDirt(MoDirt) =='mo':
            threshedImage = threshedImage.astype(np.bool)*blur.astype(np.bool)
            threshedImage =threshedImage.astype(np.uint8)*255
        return threshedImage
    elif type(Mask)==np.ndarray:
        raise Exception("Mask and image have different dimensions!")
    else:
        # If no mask, still
        if fun.checkMoDirt(MoDirt) =='dirt':
            # Invert the image so it comes out as white dirt on black
            threshedImage = np.bitwise_not(threshedImage)
            # Dirt particle counting requires an 8-bit image with white as 255
            threshedImage = threshedImage.astype(np.uint8)
            threshedImage[threshedImage!=0]=255
            return threshedImage
        elif fun.checkMoDirt(MoDirt) =='mo':
            return threshedImage

################################################################################

################################################################################

def bigAmountDirt(img):
    """
    OBSLETE. TRANSITIONED TO measure.calcDirt.
    Calculates the number of dirt particles and the area of the foil covered by dirt
    Takes a black and white image (black dirt on white foil)
    Returns the number of dirt particles, the area of dirt particles, and a labelled image
    inverse the colors since mahotas.label only works on white on black
    """
    #inv = cv2.bitwise_not(img)
    inv=img.copy()
    #invDirt = cv2.bitwise_not(isoDirt(img,profile))
    labeledFoil,numDirt = mh.label(inv)
    # Don't count the foil in numDirt:
    
    #Calculate the area of the dirt using Findcontours
    sizes = mh.labeled.labeled_size(labeledFoil)
    # Sort sizes of particles by size in descending order:
    sizes = np.sort(sizes)[::-1]
    # Eliminate the foil from the "sizes" array:
    sizes = sizes[2:]
    # Total area of dirt is equal to the sum of the sizes.
    areaDirt = sum(sizes)
    return (numDirt, areaDirt)

################################################################################

################################################################################

def makeDirtMap(panImage, maskImage, MoDirt='dirt'):
    """
    Test function
    """
    # Create the threshholded image
    subImage = fun.loadImg(panImage,0)
    subMask = fun.loadImg(maskImage,0)
    proc = images.bigPostPreProc(subImage)
    poster = images.bigPosterfy(proc)
    threshed = bigRegionalThresh(subImage,poster,Mask=subMask,\
    p=8,d=28,m=55,pt=60,hE=90,MoDirt=MoDirt)
    threshed = threshed.astype(np.uint8)*255
    return threshed, subMask

################################################################################

################################################################################

def countDirt(threshed, subMask,res=4,verbose=True):
    """
    Test function
    """
    # Get amount of dirt and write output to csv file:
    numDirt, PixDirt = bigAmountDirt(threshed)
    if verbose: print "NumDirt: "+str(numDirt)
    if verbose: print "PixDirt: " + str(PixDirt)
    AreaDirt = round(float(PixDirt)*res*10**-6, 6)
    PixFoil = np.sum(subMask.astype(np.bool_))
    AreaFoil = round(float(PixFoil)*res*10**-6, 6)
    if AreaFoil == 0:
        PercDirt = np.nan
    else:
        PercDirt = round(float(AreaDirt)/float(AreaFoil)*100,2)
    return PercDirt

################################################################################

################################################################################
# Also can be found in sandbox.foldertools. Will probably consolidate that later.
def FILonlySubimages(contents, limitToType=0):
    '''
    Takes a list of folder contents and returns only the items that are subImages
    Kwargs:
         - limitToType - if limitToType is a string, return only items ending with
                         that string, i.e. limitToType = '.tif' filters all but .tif
                         files.
                       - limitToType can also be a list, set, or tuple of strings
                       - if limitToType is anything else, it will not do anything
    '''
    contents[:] = [img for img in contents if img.startswith("sub_")]
    if type(limitToType) in (list, tuple, set):
        contents[:] = [img for img in contents if os.path.splitext(img) in limitToType]
    elif type(limitToType)==str:
        contents[:] = [img for img in contents if img.endswith(limitToType)]
    return contents

    
"""
EXTRA CODE:
       
def test(subfolder):
    contents = os.listdir(subfolder)
    data = open("data.csv","w+b")
    writer = csv.writer(data)
    writer.writerow(["Foilnum","Count","Area"])
    for sub in contents:
        image = fun.loadImg(subfolder+"/"+sub,0)
        numDirt, areaDirt = bigAmountDirt(image)
        print sub + " Count: "+str(numDirt)+" Area: "+str(areaDirt)
        writer.writerow([sub, numDirt, areaDirt])
        
"""