import cv2
import os
import numpy as np
import GenSIP.functions as fun
import GenSIP.measure as meas
import GenSIP.gencsv as gencsv


from GenSIP.cleantests.moly import Monalysis
from GenSIP.cleantests.dirt import dirtnalysis
from GenSIP.bigscans.bigfoils import ImgAnalysis



def analyzeImgOrFolder (path, res, **kwargs):
    """
    Receives a path to an image or a folder of images and saves the output csv
    file and pictures to a folder 
    Example:
        > analyzeImgOrFolder ("InputPicts/Folder/",16,MoDirt='mo',method='cleantests')
    Will run the cleantest 
    """
    method=kwargs.get('method', "cleantests")
    MoDirt=kwargs.get('MoDirt', 'Mo')
    Mask=kwargs.get('Mask', 0)
    genPoster=kwargs.get('genPoster', False)
    MoDirt = fun.checkMoDirt(MoDirt)
    
    filetypes = ['.tif', '.jpg', '.jpeg','.tiff']
    
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
        
    outFolder = "Output/Output_"+name
    
    if genPoster: posterFolder = outFolder+'/PosterMaps/'
    
    if MoDirt == 'mo':
        mapFolder = os.path.join(outFolder,'PtMaps/')
    else:
        mapFolder = os.path.join(outFolder,'DirtMaps/')
        
    if not os.path.exists(mapFolder): os.makedir(mapFolder)
    if not os.path.exists(mapFolder): os.makedir(mapFolder)
    if genPoster and not os.path.exists(posterFolder): os.makedir(posterFolder)
    
    Data = {}
    if os.path.isdir(path):
        images = [f for f in os.listdir(path) if os.path.splitext(f)[1] in filetypes]
        imgPaths = [os.path.join(path,f) for f in images]
        
        for i in range(len(images)):
            statsDict, picts = analyzeImage(imgPaths[i], res, 
                                            method=method, MoDirt=MoDirt, 
                                            Mask=Mask)
            imgName = os.path.splitext(images[i])[0]
            Data[imgName] = statsDict
            (threshed,
             poster) = picts
            # Create the output images
            cv2.imwrite(mapFolder+name+'.png',
                        threshed, [cv2.cv.CV_IMWRITE_PNG_COMPRESSION,6])
            if genPoster:
                cv2.imwrite(posterFolder+name+'.png',
                            poster, [cv2.cv.CV_IMWRITE_PNG_COMPRESSION,6])
    else:
        statsDict, picts = analyzeImage(path, res, 
                                        method=method, MoDirt=MoDirt, 
                                        Mask=Mask)
        Data[name] = statsDict
        (threshed,
            poster) = picts
        # Create the output images
        cv2.imwrite(mapFolder+name+'.png',
                    threshed, [cv2.cv.CV_IMWRITE_PNG_COMPRESSION,6])
        if genPoster:
            cv2.imwrite(posterFolder+name+'.png',
                        poster, [cv2.cv.CV_IMWRITE_PNG_COMPRESSION,6])
    filePath = os.path.join(outFolder,MoDirt.capitalize()+'_ouput_'+name+'.csv')
    CSV = gencsv.DataToCSV(filePath, name)   
    CSV.writeDataFromDict(Data,FirstColHead='Image')
    CSV.closeCSVFile()         
    
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
    
    if MoDirt == 'mo':
        if method.lower() in ['cleantests','smallfoils', 'cleantest']:
            (PtArea, 
            FoilArea, 
            MolyArea, 
            MolyMass, 
            threshed) = Monalysis(img, res,verbose=verbose)
            
            PercPt = 100*PtArea/FoilArea
            poster = fun.makePoster(img)
            
        elif method.lower() in ['bigfoils','big','bigscans','no border']:
            stats, picts = ImgAnalysis(img, mask, res, MoDirt=MoDirt,returnSizes=False)
            (PtArea,
            FoilArea,
            PercPt) = stats
            MolyArea = FoilArea-PtArea
            MolyMass = MolyArea*.3*10.2 #moly mass in micrograms
            (threshed, poster) = picts
            
        retData = {'Pt Area (mm^2)':PtArea,
                    'Foil Area (mm^2)':FoilArea,
                    'Moly Area (mm^2)':MolyArea,
                    'Mass Molybdenum (micrograms)':MolyMass,
                    '% Exposed Pt':PercPt}
                    
    elif MoDirt == 'dirt':
        if method.lower() in ['cleantests','smallfoils', 'cleantest']:
            (DirtNum,
             DirtArea,
             threshed,
             DirtSizes) = dirtnalysis (img, res, MaskEdges=True, retSizes=True)
            DirtSizes = DirtSizes*res
            MeanSize = DirtSizes.mean()
            # Number of Particles with a diameter approximately over 100 microns,
            # Corresponding to an area of ~7854 square microns. 
            numOver100 = DirtSizes[DirtSizes>7850].size
            poster = fun.makePoster(img)
            
        elif method.lower() in ['bigfoils','big','bigscans','no border']:
            stats, picts = ImgAnalysis(img, mask, res, MoDirt=MoDirt,returnSizes=True)
            (DirtNum,
             DirtArea,
             AreaFoil,
             Perc,
             DirtSizes) = stats
             
            DirtSizes = DirtSizes*res
            MeanSize = DirtSizes.mean()
            # Number of Particles with a diameter approximately over 100 microns,
            # Corresponding to an area of ~7854 square microns. 
            numOver100 = DirtSizes[DirtSizes>7850].size
            MaxSize = DirtSizes.max()
            
            (threshed, poster) = picts
            
        retData = {'Dirt Count':DirtNum,
                    'Dirt Area (mm^2)':DirtArea,
                    'Mean Particle Area (micron^2)':MeanSize,
                    'Max Particle Area (micron^2)':MaxSize,
                    '# Dirt Particles over 100micron diameter':numOver100}
    
    retPicts = (threshed,poster)
        
    return retData, retPicts
      
                             
                           
"""
            
    if method.lower() in ['cleantests','smallfoils', 'cleantest']:
        if MoDirt == 'mo':
            (PtArea, 
            FoilArea, 
            MolyArea, 
            MolyMass, 
            threshed) = Monalysis(img, res,verbose=verbose)
            
            PercPt = 100*PtArea/FoilArea
            
            retData = {'Pt Area (mm^2)':PtArea,
                       'Foil Area (mm^2)':FoilArea,
                       'Moly Area (mm^2)':MolyArea,
                       'Mass Molybdenum (micrograms)':MolyMass,
                       '% Exposed Pt':PercPt}
        elif MoDirt=='dirt':
            (DirtNum,
             DirtArea,
             threshed,
             DirtSizes) = dirtnalysis (img, res, MaskEdges=True, retSizes=True)
            DirtSizes = DirtSizes*res
            MeanSize = DirtSizes.mean()
            # Number of Particles with a diameter approximately over 100 microns,
            # Corresponding to an area of ~7854 square microns. 
            numOver100 = DirtSizes[DirtSizes>7850].size
            MaxSize = DirtSizes.max()
            retData = {'Dirt Count':DirtNum,
                       'Dirt Area (mm^2)':DirtArea,
                       'Mean Particle Area (micron^2)':MeanSize,
                       'Max Particle Area (micron^2)':MaxSize,
                       '# Dirt Particles over 100micron diameter':numOver100}
        poster = fun.makePoster(img)
    elif method.lower() in ['bigfoils','big','bigscans','no border']:
        if MoDirt == 'mo':
            stats, picts = ImgAnalysis(img, mask, res, MoDirt=MoDirt,returnSizes=False)
            (PtArea,
            FoilArea,
            PercPt) = stats
            MolyArea = FoilArea-PtArea
            MolyMass = MolyArea*.3*10.2 #moly mass in micrograms
            retData = {'Pt Area (mm^2)':PtArea,
                       'Foil Area (mm^2)':FoilArea,
                       'Moly Area (mm^2)':MolyArea,
                       'Mass Molybdenum (micrograms)':MolyMass,
                       '% Exposed Pt':PercPt}
        elif MoDirt == 'dirt':
            (DirtNum,
             DirtArea,
             AreaFoil,
             Perc,
             DirtSizes) = ImgAnalysis(img, mask, res, MoDirt=MoDirt,returnSizes=True)
             
            DirtSizes = DirtSizes*res
            MeanSize = DirtSizes.mean()
            # Number of Particles with a diameter approximately over 100 microns,
            # Corresponding to an area of ~7854 square microns. 
            numOver100 = DirtSizes[DirtSizes>7850].size
            MaxSize = DirtSizes.max()
            retData = {'Dirt Count':DirtNum,
                       'Dirt Area (mm^2)':DirtArea,
                       'Mean Particle Area (micron^2)':MeanSize,
                       'Max Particle Area (micron^2)':MaxSize,
                       '# Dirt Particles over 100micron diameter':numOver100}
"""
        