# -*- coding: utf-8 -*-
import cv2
import os
import numpy as np
import GenSIP.functions as fun
import GenSIP.measure as meas
import GenSIP.gencsv as gencsv
import GenSIP.histomethod.display as dis

from GenSIP.cleantests.moly import Monalysis
from GenSIP.cleantests.dirt import dirtnalysis
from GenSIP.bigscans.bigfoils import ImgAnalysis
from GenSIP.histomethod.mainanalysis import analyzeByHisto

################################################################################

################################################################################

def analyzImgOrFolder (path, res, method="cleantests", **kwargs):
    """
    Analyze Image or Folder (analyzImgOrFolder)
    
    Receives a path to an image or a folder of images and saves the output csv
    file and pictures to a folder in the Output/ directory. 
    
    ARGUMENTS:
        - path - path to the folder or image to be analyzed
        - res - resolution of the image, in microns^2/pixel area
        
    KEY-WORD ARGUMENTS:
        - method =  "cleantests" - the method to use to analyze the image.
                    Can be one of the following:
                        cleantests - method used by cleantests.dirt and moly
                        bigfoils - method used by bigscans.bigfoils
                        histogram - new method developed that automatically assigns
                            threshold values based on the histogram of the image.
                            Currently analyzImgOrFolder is the most refined means
                            of using the histogram module. 
                        standards - this method works by accessing the manually 
                            thresholded images stored in standards/all_dirt/ and 
                            standards/all_plat/ inside the current working 
                            directory, and calculates the platinum or dirt data
                            directly from those thresholded images. The path variable 
                            must be a folder of the standard images. 
        - MoDirt = 'Mo' - option to do molybdenum or dirt analysis
        - Mask = 0 - option to include a path string to either a single mask (in 
                    the case that the path variable links to a single image) or 
                    a folder of masks with the same name as the image they corre-
                    spond to. i.e. the image is "Images/sub_img_004_003.png" and the 
                    mask is at "masks/sub_img_004_003.png"
        - genPoster = False - Option to generate a folder of poster images in the
                    the output folder.
        - autoMaskEdges = False - allows the bigfoil and histogram methods to work
                    with individual foils like the cleantest method by automatically 
                    masking off the dark background around the foil.If this option
                    is set to True, then any input for the Mask variable is over
                    -ridden.
        - stdDir = "standards/" - sets the directory containing the standards in
                    the "all_stds/", "all_plat/", "all_dirt/" and "all_masks" directories.
        - verbose = False - makes the function verbose.
        - compToStds = False - Currently not used for anything. Would like to 
                    turn this option into a means of directly comparing the 
                    results of one method to the standard results

    EXAMPLE:
        > analyzImgOrFolder ("InputPicts/Folder/",16,MoDirt='mo',method='cleantests')
    Will run the cleantest analysis on all images in "InputPicts/Folder/" and 
    will save all the output to a folder "Output/Output_Folder_cleantests/". 
    The contents of this folder are:
        - PtMaps/ - a folder of the thresholded images of exposed platinum
        - Mo_output_Folder.csv - a CSV file that contains the molybdenum analysis data

    """
    MoDirt=kwargs.get('MoDirt', 'Mo')
    Mask=kwargs.get('Mask', 0)
    genPoster=kwargs.get('genPoster', False)
    # compareToStds = kwargs.get('compToStds',False)
    verbose = kwargs.get('verbose',False)
    autoMask = kwargs.get('autoMaskEdges',False)
    stdDir = kwargs.get('stdDir', 'standards/')
    
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
    outFolder = "Output/Output_"+name+'_'+method
    
    if genPoster: posterFolder = outFolder+'/PosterMaps/'
    
    if MoDirt == 'mo':
        mapFolder = os.path.join(outFolder,'PtMaps/')
    else:
        mapFolder = os.path.join(outFolder,'DirtMaps/')
        
    if not os.path.exists(mapFolder): os.makedirs(mapFolder)
    if not os.path.exists(mapFolder): os.makedirs(mapFolder)
    if genPoster and not os.path.exists(posterFolder): os.makedirs(posterFolder)
    
    # Verbose Feedback:
    if verbose: 
        if os.path.isdir(path): foo = "directory"
        else: foo = "image"
        print "––––––––––––NEXUS analyzImgOrFolder––––––––––––––"
        print "Now Running " + MoDirt + " testing using the " + method + " method"
        print "On "+foo+": '"+name+"':"
        print "  Resolution set to: " + str(res)
        print "  Mask set to: "+ str(Mask)
        print "  autoMaskEdges set to: "+str(autoMask)
        if method in ['standards','standard','std','stds']:
            print "  Standard Directory: " + str(stdDir)
        print "- - - - - - - - - - - - - - - - - - - - - - - - - - - -"
    #--------------#
    # Debug switch #
    debug = 0      #
    #--------------#
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
            if debug: 
                for m in masks: print m
            # Create paths to those images
            maskPaths = [os.path.join(Mask,m) for m in masks]
            if debug: 
                for i in range(len(maskPaths)): 
                    print maskPaths[i] + str(fun.loadImg(maskPaths[i]).shape)
                    print imgPaths[i] + str(fun.loadImg(imgPaths[i]).shape)
            # I am assuming the mask name will be the same as the corresponding 
            # name in the image folder, so when both are sorted, they should match. 
            maskPaths.sort() 
            
        else:
            maskPaths = [0 for f in imgPaths]
        
        for i in range(len(images)):
            # Make the mask image from the mask path
            if Mask!=0: mask = fun.loadImg(maskPaths[i])
            else: mask=0
            imgName = os.path.splitext(images[i])[0]
            # run analysis on the image
            statsDict, picts = analyzeImage(imgPaths[i], res, 
                                            method=method, MoDirt=MoDirt, 
                                            Mask=mask,autoMaskEdges=autoMask,
                                            stdDir=stdDir, 
                                            name=imgName,
                                            outDir=outFolder,
                                            verbose=verbose)
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
        if Mask!=0 and type(Mask)==str: mask = fun.loadImg(Mask)
        elif type(Mask)==np.ndarray: mask = Mask.copy()
        else: mask=0
        statsDict, picts = analyzeImage(path, res, 
                                        method=method, MoDirt=MoDirt, 
                                        Mask=mask,autoMaskEdges=autoMask,
                                        stdDir=stdDir,
                                        name=name, 
                                        outDir=outFolder,
                                        verbose=verbose)
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
    CSV = gencsv.DataToCSV(filePath, name+': '+method+' method of analysis')   
    CSV.writeDataFromDict(Data,FirstColHead='Image')
    CSV.closeCSVFile() 
            
################################################################################

################################################################################

def analyzeImage(path, res, method='cleantests', MoDirt='mo', Mask=0, name='',
                 autoMaskEdges=False, stdDir='standards/', outDir="Output/",verbose=False):
    """
    Given the path, runs analysis on a single image using one of the methods in
    GenSIP specified by the 'method' kwarg (currently: cleantests or bigfoils). 
    Returns a Data Dictionary and the thresholded image and poster.
    """
    
    img = fun.loadImg(path)
    MoDirt = fun.checkMoDirt(MoDirt)
    
    if type(Mask) in [int,float] and Mask==0:
        mask = np.ones(img.shape)
    elif type(Mask)==np.ndarray and Mask.shape == img.shape:
        mask = Mask.copy()
    else:
        if type(Mask)==np.ndarray:
            MaskDims = str(Mask.shape)
            ImgDims = str(img.shape)
        else:
            MaskDims = '--'
            ImgDims = '--'
        raise Exception ("""There is a problem with the mask: \n
                            Mask type: {0} \n
                            Mask shape: {1} \n
                            Image shape: {2}""".format(str(type(Mask)),MaskDims,ImgDims))
    # Uses my OLD maskEdges function to mask off the dark area around a foil if 
    # specified. 
    if autoMaskEdges:
        maskedImg, mask = fun.maskEdge(img)

    # Verbose Feedback
    if verbose:
        print "Now analyzing Image at:"
        print str(path)
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
        elif method.lower() in ['newmethod', 'histograms', 'histo', 'histogram', 'hist']:
            stats, picts = analyzeByHisto (img, res, 
                                           Mask=mask, verbose=verbose,
                                           MoDirt=MoDirt, returnPoster=True,
                                           returnData=True,returnSizes=False)
            (PtArea,
            PercPt,
            FoilArea,
            Data) = stats
            
            dis.saveHist(Data, outDir, name=name)
            
            MolyArea = FoilArea-PtArea
            MolyMass = MolyArea*.3*10.2 #moly mass in micrograms
            
            (threshed, poster) = picts
            
        # STANDARD ANALYSIS –––––––––––––––––––––––––––––––––––––––––––
        elif method.lower() in ['standards','standard','std','stds']:
            poster = fun.posterfy(img)
            imgName = os.path.splitext(os.path.split(path)[1])[0]
            PtMapPath = os.path.join(stdDir, 'all_plat/')+imgName+'.png'
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
        retData = {'Pt Area (mm^2)':ifFloatRound(PtArea,4),
                    'Foil Area (mm^2)':ifFloatRound(FoilArea,2),
                    'Moly Area (mm^2)':ifFloatRound(MolyArea,3),
                    'Mass Molybdenum (micrograms)':ifFloatRound(MolyMass,3),
                    '% Exposed Pt':ifFloatRound(PercPt,3)}
                    
    # DIRT ANALYSIS ============================================================
    elif MoDirt == 'dirt':
        
        # Method used by cleantests  ––––––––––––––––––––––––––––––––––––
        if method.lower() in ['cleantests','smallfoils', 'cleantest','ct','foils']:
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
                                           returnData=True,returnSizes=True)
            (DirtNum,
             DirtArea,
             DirtSizes,
             AreaFoil,
             Data) = stats
            
            dis.saveHist(Data, outDir, name=name)
            
            (threshed, poster) = picts
            
        # STANDARD ANALYSIS –––––––––––––––––––––––––––––––––––––––––––
        elif method.lower() in ['standards','standard','std','stds']:
            poster = fun.posterfy(img)
            imgName = os.path.splitext(os.path.split(path)[1])[0]
            DirtMapPath = os.path.join(stdDir, 'all_dirt/')+imgName+'.png'
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

        retData = {'Dirt Count':ifFloatRound(DirtNum, 2),
                    'Dirt Area (mm^2)':ifFloatRound(DirtArea, 5),
                    'Mean Particle Area (micron^2)':ifFloatRound(MeanSize,1),
                    'Max Particle Area (micron^2)':ifFloatRound(MaxSize,1),
                    '% Dirt Particles over 100micron diameter':ifFloatRound(percOver100,3)}
    if verbose:
        print "- - - - - - - - - - - - - - - - - - - - - - - - - - - -"
    
    # Return results
    retPicts = (threshed,poster)
        
    return retData, retPicts
      
################################################################################

################################################################################

def ifFloatRound(num, rnd):
    """ 
    If num is a float, returns num rounded by rnd digits.
    Otherwise just return num. 
    Prevents errors by rounding strings
    """
    if type(num) == float: return round(num,rnd)
    
    else: return num

################################################################################

################################################################################

def blankDataDict(MoDirt='mo'):
    """
    Returns a data dictionary of the format used in analyzImgOrFolder where all 
    of the data points are replaced with "'--" as a filler."
    """
    MoDirt = fun.checkMoDirt(MoDirt)
    if MoDirt=='mo':
        """Gives a blank moly data dictionary"""
        retData = {'Pt Area (mm^2)':"'--",
                    'Foil Area (mm^2)':"'--",
                    'Moly Area (mm^2)':"'--",
                    'Mass Molybdenum (micrograms)':"'--",
                    '% Exposed Pt':"'--"}
    elif MoDirt=='dirt':
        """Gives a blank dirt data dictionary"""
        retData = {'Dirt Count':"'--",
                    'Dirt Area (mm^2)':"'--",
                    'Mean Particle Area (micron^2)':"'--",
                    'Max Particle Area (micron^2)':"'--",
                    '% Dirt Particles over 100micron diameter':"'--"}
    return retData

################################################################################

################################################################################

def blankImg(dims=(100,100)):
    """Creates a white image given the dimensions"""
    return np.ones(dims)*255

################################################################################

################################################################################
# 
def compareToStandards(method, res, STDpath='standards/all_stds'):
    """Creates a csv file of the standard data"""
    dirtFolder = os.path.join(STDpath,'all_dirt')
    platFolder = os.path.join(STDpath,'all_plat')
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
'''
#TESTING:
x = raw_input("Are you going to run the full test? (y/n) ")
if x=='n':# on-off switch for running this test.
    pass
elif x=='y':
    if __name__=='__main__':
        import shutil
        path = "standards/all_stds"
        methods = ['cleantests','bigfoils', 'histogram','standard']
        masks = [0,"standards/all_masks"]
        modirt = ['Mo','Dirt']
        for meth in methods:
            for m in masks:
                for md in modirt:
                    analyzImgOrFolder(path, 16, method=meth,
                                    MoDirt=md,Mask=m,
                                    genPoster=True,verbose=True)
                    if type(m)==str:
                        M = 'withMasks'
                    else:
                        M = '0Masks'
                    direct = "_".join(["Output/Tests/master",meth, M, md])
                    os.makedirs(direct)
                    shutil.move("Output/Output_all_stds_"+meth,direct)
                    
        folder = "Output/Tests/"
        if not os.path.exists("Output/Tests/all_stds"):
            os.makedirs("Output/Tests/all_stds")
        for direct in [d for d in os.listdir(folder) if d.startswith('master')]:
            meth = direct.split('_')[1]
            path = os.path.join(folder,direct,"Output_all_stds_"+meth)
            CSVname = [c for c in os.listdir(path) if c.endswith('.csv')][0]
            CSVnewpath = os.path.join("Output/Tests/all_stds",direct+'.csv')
            CSVoldpath = [os.path.join(path,CSVname)][0]
            print "old path: " + CSVoldpath
            print "new path: " + CSVnewpath
            print "name: " + CSVname
            if os.path.exists(CSVoldpath):
                os.rename(CSVoldpath,CSVnewpath)

'''