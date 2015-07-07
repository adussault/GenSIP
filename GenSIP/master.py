# -*- coding: utf-8 -*-
import cv2
import os
import numpy as np
import GenSIP.functions as fun
import GenSIP.measure as meas
import GenSIP.gencsv as gencsv

from GenSIP.cleantests.moly import Monalysis
from GenSIP.cleantests.dirt import dirtnalysis
from GenSIP.bigscans.bigfoils import ImgAnalysis
from GenSIP.histomethod.mainanalysis import analyzeByHisto

################################################################################

################################################################################

def analyzeImgOrFolder (path, res, **kwargs):
    """
    Receives a path to an image or a folder of images and saves the output csv
    file and pictures to a folder 
    Example:
        > analyzeImgOrFolder ("InputPicts/Folder/",16,MoDirt='mo',method='cleantests')
    Will run the cleantest 
    
        method=kwargs.get('method', "cleantests")
    MoDirt=kwargs.get('MoDirt', 'Mo')
    Mask=kwargs.get('Mask', 0)
    genPoster=kwargs.get('genPoster', False)
    compareToStds = kwargs.get('compToStds',False)
    verbose = kwargs.get('verbose',False)

    """
    method=kwargs.get('method', "cleantests")
    MoDirt=kwargs.get('MoDirt', 'Mo')
    Mask=kwargs.get('Mask', 0)
    genPoster=kwargs.get('genPoster', False)
    compareToStds = kwargs.get('compToStds',False)
    verbose = kwargs.get('verbose',False)
    
    # Standardize MoDirt to 'mo' or 'dirt' using checkMoDirt
    MoDirt = fun.checkMoDirt(MoDirt)
    
    filetypes = ['.tif', '.jpg', '.jpeg','.tiff']
    
    # Standardize the path string, and extract the name of the folder or image 
    # file depending on whether the path is the path to a directory or image. 
    if os.path.isdir(path):
        if path.endswith('/'):
            path = path[:-1]
        name = os.path.split(path)[1]
    elif os.path.splitext(path)[1] in filetypes:
        name = os.path.splitext(os.path.split(path)[1])[0]
    elif type(path)!=str: 
        raise Exception("Path must be a string: %s" % str(path))
    else: 
        raise Exception("Invalid path name: %s" % path)
        
    # Generate output folders
    outFolder = "Output/Output_"+name
    
    if genPoster: posterFolder = outFolder+'/PosterMaps/'
    
    if MoDirt == 'mo':
        mapFolder = os.path.join(outFolder,'PtMaps/')
    else:
        mapFolder = os.path.join(outFolder,'DirtMaps/')
        
    if not os.path.exists(mapFolder): os.makedirs(mapFolder)
    if not os.path.exists(mapFolder): os.makedirs(mapFolder)
    if genPoster and not os.path.exists(posterFolder): os.makedirs(posterFolder)
    
    """Create Data Dictionary"""
    # Iterate through the images within a folder if the path is to a directory, 
    # and run analyzeImg on each of image, then write the results to the Data 
    # Dictionary. 
    Data = {}
    
    # OPERATE ON FOLDER OF IMAGES ==============================================
    if os.path.isdir(path):
        
        # Get list of images in directory
        images = [f for f in os.listdir(path) if os.path.splitext(f)[1] in filetypes]
        # Create paths to those images
        imgPaths = [os.path.join(path,f) for f in images]
        imgPaths.sort()
        if Mask!=0:
            assert type(Mask)==str, """
                                    'Mask' kwarg must be a path to a directory
                                    if the 'path' variable is a path to a directory."
                                    """
            assert os.path.isdir(Mask), """
                                        'Mask' kwarg must be a path to a directory
                                        if the 'path' variable is a path to a directory.
                                        """
            # Get list of images in directory
            masks  = [m for m in os.listdir(Mask) if os.path.splitext(m)[1] in filetypes]
            # Create paths to those images
            maskPaths = [os.path.join(Mask,m) for m in masks]
            # I am assuming the mask name will be the same as the corresponding 
            # name in the image folder, so when both are sorted, they should match. 
            maskPaths.sort() 
            
        else:
            maskPaths = [0 for f in imgPaths]
        
        for i in range(len(images)):
            # Make the mask image from the mask path
            if Mask!=0: mask = fun.loadImg(maskPaths[i])
            else: mask=0
            # run analysis on the image
            statsDict, picts = analyzeImage(imgPaths[i], res, 
                                            method=method, MoDirt=MoDirt, 
                                            Mask=mask)
            imgName = os.path.splitext(images[i])[0]
            # Assign to Data Dictionary
            Data[imgName] = statsDict
            (threshed,
             poster) = picts
            threshed = threshed.astype(np.uint8)
            threshed[threshed!=0]=255
            poster = poster.astype(np.uint8)
            
            # Create the output images
            cv2.imwrite(mapFolder+imgName+'.png',
                        threshed, [cv2.cv.CV_IMWRITE_PNG_COMPRESSION,6])
            if genPoster:
                cv2.imwrite(posterFolder+imgName+'.png',
                            poster, [cv2.cv.CV_IMWRITE_PNG_COMPRESSION,6])
                            
    # OPERATE ON A SINGLE IMAGE ================================================
    else:
        # run analysis on the image
        statsDict, picts = analyzeImage(path, res, 
                                        method=method, MoDirt=MoDirt, 
                                        Mask=Mask)
        Data[name] = statsDict
        (threshed,
         poster) = picts
        threshed = threshed.astype(np.uint8)
        threshed[threshed!=0]=255
        poster = poster.astype(np.uint8)
        poster[poster!=0]=255
        # Create the output images
        cv2.imwrite(mapFolder+name+'.png',
                    threshed, [cv2.cv.CV_IMWRITE_PNG_COMPRESSION,6])
        if genPoster:
            cv2.imwrite(posterFolder+name+'.png',
                        poster, [cv2.cv.CV_IMWRITE_PNG_COMPRESSION,6])
                        
    """Write the output to a CSV file"""
    filePath = os.path.join(outFolder,MoDirt.capitalize()+'_ouput_'+name+'.csv')
    CSV = gencsv.DataToCSV(filePath, name)   
    CSV.writeDataFromDict(Data,FirstColHead='Image')
    CSV.closeCSVFile() 
            
################################################################################

################################################################################

def analyzeImage(path, res, method='cleantests', MoDirt='mo', Mask=0, verbose=False):
    """
    Given the path, runs analysis on a single image using one of the methods in
    GenSIP specified by the 'method' kwarg (currently: cleantests or bigfoils). 
    Returns a Data Dictionary and the thresholded image and poster.
    """
    img = fun.loadImg(path)
    MoDirt = fun.checkMoDirt(MoDirt)
    
    if Mask==0:
        mask = np.ones(img.shape)
    elif type(Mask)==np.ndarray and Mask.shape == img.shape:
        mask = Mask.copy()
    else:
        raise Exception
    retData = {}
    
    # MOLYBDENUM ANALYSIS ======================================================
    if MoDirt == 'mo':
        # Method used by cleantests  ––––––––––––––––––––––––––––––––––––
        if method.lower() in ['cleantests','smallfoils', 'cleantest']:
            (PtArea, 
            FoilArea, 
            MolyArea, 
            MolyMass, 
            threshed) = Monalysis(img, res,verbose=verbose)
            
            PercPt = 100*PtArea/FoilArea
            poster = fun.makePoster(img)
        
        # Method used by bigfoils  –––––––––––––––––––––––––––––––––––––
        elif method.lower() in ['bigfoils','big','bigscans','no border']:
            stats, picts = ImgAnalysis(img, mask, res, MoDirt=MoDirt,returnSizes=False)
            (PtArea,
            FoilArea,
            PercPt) = stats
            MolyArea = FoilArea-PtArea
            MolyMass = MolyArea*.3*10.2 #moly mass in micrograms
            (threshed, poster) = picts
            
        # Method Used by Histogram Analysis (newmethod) ––––––––––––––––––––
        elif method.lower() in ['newmethod', 'histograms', 'histo', 'histogram']:
            stats, picts = analyzeByHisto (img, res, 
                                           Mask=mask, verbose=verbose,
                                           MoDirt=MoDirt, returnPoster=True,
                                           returnData=False,returnSizes=False)
            (PtArea,
            PercPt,
            FoilArea) = stats
            
            MolyArea = FoilArea-PtArea
            MolyMass = MolyArea*.3*10.2 #moly mass in micrograms
            
            (threshed, poster) = picts
            
        # STANDARD ANALYSIS –––––––––––––––––––––––––––––––––––––––––––
        elif method.lower() in ['standards','standard','std','stds']:
            poster = fun.posterfy(img)
            imgName = os.path.splitext(os.path.split(path)[1])[0]
            PtMapPath = 'standards/all_plat/'+imgName+'.png'
            if os.path.exists(PtMapPath):
                threshed = fun.loadImg(PtMapPath)
                PtArea = meas.calcExposedPt(threshed, res, getAreaInSquaremm=True)
                PixFoil = np.sum(mask.astype(np.bool_))
                FoilArea = round(PixFoil*res*10**-6, 4)
                MolyArea = FoilArea-PtArea
                MolyMass = MolyArea*.3*10.2 #moly mass in micrograms
                if FoilArea == 0:
                    PercPt = 0
                else:
                    PercPt = round(float(PtArea)/float(FoilArea)*100,2)
            else:
                print "Not a standard: " + imgName
                print "  File path does not Exist: " + PtMapPath
                retData = blankDataDict(MoDirt)
                threshed = blankImg(img.shape)
                return retData, (threshed, poster)
                
        # UNMATCHED METHOD ––––––––––––––––––––––––––––––––––––––––––––
        else:
            raise Exception("""The specified method is not available: {0} \n
                               Method should be one of the following: \n
                               'cleantests','bigfoils','histogram','standard'.
                               """.format(str(method)))
                            
        # Prepare Return Data Dictionary ---------------------------------------
        retData = {'Pt Area (mm^2)':round(PtArea,4),
                    'Foil Area (mm^2)':round(FoilArea,2),
                    'Moly Area (mm^2)':round(MolyArea,3),
                    'Mass Molybdenum (micrograms)':round(MolyMass,3),
                    '% Exposed Pt':round(PercPt,3)}
                    
    # DIRT ANALYSIS ============================================================
    elif MoDirt == 'dirt':
        
        # Method used by cleantests  ––––––––––––––––––––––––––––––––––––
        if method.lower() in ['cleantests','smallfoils', 'cleantest']:
            (DirtNum,
             DirtArea,
             threshed,
             DirtSizes) = dirtnalysis (img, res, MaskEdges=True, retSizes=True)
                         
            poster = fun.makePoster(img)
            
        # Method used by bigfoils  –––––––––––––––––––––––––––––––––––––
        elif method.lower() in ['bigfoils','big','bigscans','no border']:
            stats, picts = ImgAnalysis(img, mask, res, 
                                       MoDirt=MoDirt,returnSizes=True)
            (DirtNum,
             DirtArea,
             AreaFoil,
             Perc,
             DirtSizes) = stats
                                     
            (threshed, poster) = picts
            
        # Method Used by Histogram Analysis (newmethod) ––––––––––––––––––––
        elif method.lower() in ['newmethod', 'histograms', 'histo', 'histogram']:
            stats, picts = analyzeByHisto (img, res, 
                                           Mask=mask, verbose=verbose,
                                           MoDirt=MoDirt, returnPoster=True,
                                           returnData=False,returnSizes=True)
            (DirtNum,
             DirtArea,
             DirtSizes,
             AreaFoil) = stats
            
            
            (threshed, poster) = picts
            
        # STANDARD ANALYSIS –––––––––––––––––––––––––––––––––––––––––––
        elif method.lower() in ['standards','standard','std','stds']:
            poster = fun.posterfy(img)
            imgName = os.path.splitext(os.path.split(path)[1])[0]
            DirtMapPath = 'standards/all_dirt/'+imgName+'.png'
            if os.path.exists(DirtMapPath):
                threshed = fun.loadImg(DirtMapPath)
                (DirtArea, 
                DirtNum,
                DirtSizes,
                labeled) = meas.calcDirt(threshed,
                                         res, 
                                         returnSizes=True,
                                         returnLabelled=True, 
                                         getAreaInSquaremm=True) 
            else:
                print "Not a standard: " + imgName
                print "  File path does not Exist: " + DirtMapPath
                retData = blankDataDict(MoDirt)
                threshed = blankImg(img.shape)
                return retData, (threshed, poster)
                
        # UNMATCHED METHOD ––––––––––––––––––––––––––––––––––––––––––––
        else:
            raise Exception("""The specified method is not available: {0} \n
                               Method should be one of the following: \n
                               'cleantests','bigfoils','histogram','standard'.
                               """.format(str(method)))
        
        # Prepare Return Data Dictionary ---------------------------------------
        (MeanSize, 
        MaxSize, 
        percOver100) = meas.getDirtSizeData(DirtSizes, res)
        
        retData = {'Dirt Count':DirtNum,
                    'Dirt Area (mm^2)':round(DirtArea, 5),
                    'Mean Particle Area (micron^2)':round(MeanSize,1),
                    'Max Particle Area (micron^2)':round(MaxSize,1),
                    '% Dirt Particles over 100micron diameter':round(percOver100,3)}
    
    # Return results
    retPicts = (threshed,poster)
        
    return retData, retPicts
      
################################################################################

################################################################################

def blankDataDict(MoDirt='mo'):
    MoDirt = fun.checkMoDirt(MoDirt)
    if MoDirt=='mo':
        retData = {'Dirt Count':"'--",
                'Dirt Area (mm^2)':"'--",
                'Mean Particle Area (micron^2)':"'--",
                'Max Particle Area (micron^2)':"'--",
                '% Dirt Particles over 100micron diameter':"'--"}
    elif MoDirt=='dirt':
        retData = {'Dirt Count':"'--",
                    'Dirt Area (mm^2)':"'--",
                    'Mean Particle Area (micron^2)':"'--",
                    'Max Particle Area (micron^2)':"'--",
                    '% Dirt Particles over 100micron diameter':"'--"}
    return retData

################################################################################

################################################################################

def blankImg(dims=(100,100)):
    return np.ones(dims)*255

################################################################################

################################################################################

def compareToStandards(method, res, STDpath='standards/all_stds'):
    """Creates a csv file of the standard data"""
        dirtFolder = os.path.join(path,'all_dirt')
    platFolder = os.path.join(path,'all_plat')
    dirtPaths = [os.path.join(dirtFolder,dm)
                 for dm in os.listdir(dirtFolder)
                 if dm.endswith('.png')]
    PtPaths = [os.path.join(platFolder,dm)
              for dm in os.listdir(platFolder)
              if dm.endswith('.png')]
    for path in dirtFolder:
        dirtmap = fun.loadImg(path)
        (dirtArea, 
        dirtNum,
        DirtSizes,
        labeled) = meas.calcDirt(dirtmap,
                                res, 
                                returnSizes=True,
                                returnLabelled=True, 
                                getAreaInSquaremm=True) 
        (MeanSize, 
        MaxSize, 
        percOver100) = meas.getDirtSizeData(DirtSizes, res)

################################################################################

################################################################################
                
#TESTING:
if __name__=='__main__':
    import shutil
    path = "standards/all_stds"
    methods = ['cleantests','bigfoils', 'newmethod','standard']
    masks = [0,"standards/all_masks"]
    modirt = ['Mo','Dirt']
    for meth in methods:
        for m in masks:
            for md in modirt:
                analyzeImgOrFolder(path, 16, method=meth,
                                   MoDirt=md,Maks=m,
                                   genPoster=True,verbose=True)
                if type(m)==str:
                    M = 'withMasks'
                else:
                    M = '0Masks'
                direct = "_".join(["Output/Tests/master",meth, M, md])
                os.makedirs(direct)
                shutil.move("Output/Output_all_stds",direct)
                
    folder = "Output/Tests/"
    if not os.path.exists("Output/Tests/all_stds"):
        os.makedirs("Output/Tests/all_stds")
    for direct in [d for d in os.listdir(folder) if d.startswith('master')]:
        path = os.path.join(folder,direct,"Output_all_stds")
        CSVname = [c for c in os.listdir(path) if c.endswith('.csv')][0]
        CSVnewpath = os.path.join("Output/Tests/all_stds",direct+'.csv')
        CSVoldpath = [os.path.join(path,CSVname)][0]
        print "old path: " + CSVoldpath
        print "new path: " + CSVnewpath
        print "name: " + CSVname
        if os.path.exists(CSVoldpath):
            os.rename(CSVoldpath,CSVnewpath)