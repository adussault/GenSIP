"""
Contains functions that measure the exposed platinum, dirt particles, molybdenum,
etc. based on the thresholded images produced by analysis methods.
"""
import numpy as np
import cv2
import mahotas as mh
import os

import GenSIP.functions as fun


####################################################################################

####################################################################################
             

def calcExposedPt (Ptimage, res):
    """
    Returns the area of exposed platinum in square millimeters.
    Input is the thresholded platinum image and the image resolution, 
    in square microns per pixel. 
    """
    # Convert to uint8 format
    plat = Ptimage.astype(np.uint8)
    # Make sure the image is just made of 1s and 0s
    plat[plat!=0]=1
    # Area of Pt = sum of nonzero pixels x resolution x 10^-16
    areaPt = float(plat.sum())*res*10**-6
    return areaPt

####################################################################################

####################################################################################
    
def calcDirt(img, res, returnSizes=False,returnLabelled=False,getAreaInSquaremm=False):
    """
    Calculates the number of dirt particles and the area of the foil covered by dirt
    Takes a black and white image (white dirt on black foil)
    Returns the number of dirt particles, the area of dirt particles, and a labelled image
    inverse the colors since mahotas.label only works on white on black
    """
    #inv = cv2.bitwise_not(img)
    #invDirt = cv2.bitwise_not(isoDirt(img,profile))
    labeledFoil,numDirt = mh.label(img)
    # Don't count the foil in numDirt:
    numDirt -= 1
    #Calculate the area of the dirt using Findcontours
    sizes = mh.labeled.labeled_size(labeledFoil)
    # Sort sizes of particles by size in descending order:
    sizes = np.sort(sizes)[::-1]
    # Eliminate the foil from the "sizes" array:
    sizes = sizes[2:]
    # Total area of dirt is equal to the sum of the sizes. 
    # In square microns unless otherwise specified.
    areaDirt = sum(sizes)*res
    # If specified, convert area to square mm
    if getAreaInSquaremm:
        areaDirt = float(areaDirt*10**-6)
    # Generate return tuple:
    ret = [areaDirt, numDirt]
    if returnSizes:
        ret.append(sizes)
    if returnLabelled:
        ret.append(labeledFoil)
    return tuple(ret)

####################################################################################

####################################################################################

def compareToStandards(function, res, **kwargs):
    """
    Takes a function that produces the platinum or dirt map of an image, calculates 
    the area of exposed platinum or dirt, and then compares it to manually made
    standards stored in the folder 'standards.'
    Assumes the arguments of function are (image, res, **kwargs) and the function
    returns either a black and white image of the platinum or dirt map (white on
    black), or a tuple where the first element is the platinum/dirt map and the 
    second element is the regions dictionary that contains all the histogram, 
    threshold, and morphology information for each poster region.
    
    Output format:
        The output is a dictionary with a format that varies slightly.
    
    """
    MoDirt = kwargs.get("MoDirt",'Mo')
    # Pop the values of key-word arguments that only compareToStandards uses, so
    # the remaining kwargs can be passed on to function.
    returnImages = kwargs.pop("retImages",False)
    standards = kwargs.pop("stdsDirectory",'standards/')
    usemasks = kwargs.pop("useMasks",False)
    
    directs = os.listdir(standards)
    directs = [d for d in directs if d.startswith('sub_')]
    ret = {}
    ret['MoDirt']=MoDirt
    for f in directs:
        pathtodir = standards+f
        if os.path.isdir(pathtodir) and os.path.exists(pathtodir+'/thresholds.txt'):
            try: testImg = fun.loadImg(pathtodir+'/'+f+'.tif')
            except: print pathtodir+'/'+f+'.tif'
            ret[f] = {}
            
            ### MOLYBDENUM COMPARISON________________________________________________________________
            #########################################################################################
            
            if (fun.checkMoDirt(MoDirt)=='mo') and (os.path.exists(pathtodir+'/plat.png')):
                print f
                stdPlat = fun.loadImg(pathtodir+'/plat.png')
                stdPlat[stdPlat!=0]=1
                
                if usemasks:
                    try: mask = fun.loadImg(pathtodir+'/mask.png')
                    except: mask = 0 
                    testPlat = function(testImg,res,Mask=mask,**kwargs)
                else:
                    testPlat = function(testImg,res,**kwargs)
                    
                # If the function output is a tuple, see if it is possible to 
                # extract a regions dictionary from the output and set it to regDirt.
                # Then, set testDirt to the first item in the output tuple. 
                # Otherwise, make regDict an empty dictionary and set testDirt to the 
                # output. 
                if type(testPlat)==tuple:
                    try: regDict = [d for d in testPlat if type(d)==dict][0]
                    except: regDict = {}
                    
                    testPlat = testPlat[0]

                assert type(testPlat)==np.ndarray,"Function must produce an output of type numpy.ndarray. \n"\
                +"Function output is of type {0}.".format(testPlat.__class__.__name__)
                testPlat[testPlat!=0] = 1
                testAreaPlat = calcExposedPt(testPlat,res)
                stdAreaPlat = calcExposedPt(stdPlat,res)
                fractCorrect = float(testAreaPlat)/float(stdAreaPlat)
                error = 1-fractCorrect
                
                if returnImages:
                    ret[f]['testImg'] = testPlat*255
                    ret[f]['stdImg'] = stdPlat*255

                # if regDict is the regions dictionary returned by the function,
                # incorporate the data into the return dictionary.
                if len(regDict)>0:
                    ret[f]['testData']={}
                    for reg in regDict:
                        ret[f]['testData'][reg]=regDict[reg]
                
                # Incorporate the thresholding and processing information from 
                # the thresholds.txt file that comes with the standard.
                try: threshFile = open(pathtodir+'/thresholds.txt')
                except: threshFile = 0
                if type(threshFile)==file:
                    ret[f]['stdData'] = {}
                    lines = threshFile.readlines()
                    PlatData = [line.strip().split('==') for line in lines if line.startswith(('Plat','plat'))]
                    ret[f]['stdData']['PlatProcInfo'] = [line.remove('~ ','\n') for line in lines if line.startswith('~ ')]
                    ret[f]['stdData']['PlatMaxThresh'] = PlatData[0][1]
                    ret[f]['stdData']['PlatMinThresh'] = PlatData[1][1]

                ret[f]['testAreaPt'] = testAreaPlat
                ret[f]['stdAreaPt'] = stdAreaPlat
                ret[f]['FractCorrect'] = fractCorrect
                ret[f]['FractError'] = error
               
            ### DIRT COMPARISON______________________________________________________________________
            #########################################################################################
            
            elif (fun.checkMoDirt(MoDirt)=='dirt') and (os.path.exists(pathtodir+'/dirt.png')):
                
                stdDirt = fun.loadImg(pathtodir+'/dirt.png')
                stdDirt[stdDirt!=0]=1
                if usemasks:
                    try: mask = fun.loadImg(pathtodir+'/mask.png')
                    except: mask = 0 
                    testDirt = function(testImg,res,Mask=mask,**kwargs)#STANDARD: Mask kwarg is either 0 or the mask image. 
                else:
                    testDirt = function(testImg,res,**kwargs)
                # If the function output is a tuple, see if it is possible to 
                # extract a regions dictionary from the output and set it to regDirt.
                # Then, set testDirt to the first item in the output tuple. 
                # Otherwise, make regDict an empty dictionary and set testDirt to the 
                # output. 
                try: regDict = [d for d in testDirt if type(d)==dict][0]
                except: regDict = {}
                if type(testDirt)==tuple:
                    testDirt = testDirt[0]
                
                assert type(testDirt)==np.ndarray,"Function must produce an output of type numpy.ndarray. \n"\
                +"Function output is of type {0}.".format(testDirt.__class__.__name__)
                testArea,testNum,testSizes = calcDirt(testDirt, res, returnSizes=True)
                stdinvDirt = np.bitwise_not(stdDirt)
                stdinvDirt[stdinvDirt!=0] = 1
                stdArea,stdNum,stdSizes = calcDirt(stdDirt, res, returnSizes=True)
                
                try: AreafractCorrect = float(testArea)/float(stdArea)
                except ZeroDivisionError: AreafractCorrect = np.nan
                Areaerror = 1-AreafractCorrect
                try: NumfractCorrect = float(testNum)/float(stdNum)
                except ZeroDivisionError: NumfractCorrect = np.nan
                Numerror = 1-NumfractCorrect
                
                
                ## BUILD THE RETURN DICTIONARY FOR THIS STANDARD
                if returnImages:
                    ret[f]['testImg'] = testDirt
                    ret[f]['stdImg'] = stdDirt
                    
                # if regDict is the regions dictionary returned by the function,
                # incorporate the data into the return dictionary.
                if len(regDict)>0:
                    ret[f]['testData']={}
                    for reg in regDict:
                        ret[f]['testData'][reg]=regDict[reg]
                
                # Incorporate the thresholding and processing information from 
                # the thresholds.txt file that comes with the standard.
                try: threshFile = open(pathtodir+'/thresholds.txt')
                except: threshFile = 0
                if type(threshFile)==file:
                    ret[f]['stdData'] = {}
                    lines = threshFile.readlines()
                    DirtData = [line.strip().split('==') for line in lines if line.startswith(('Dirt','dirt'))]
                    ret[f]['stdData']['DirtProcInfo'] = [line.remove('# ','\n') for line in lines if line.startswith('# ')]
                    ret[f]['stdData']['DirtMaxThresh'] = DirtData[0][1]
                    ret[f]['stdData']['DirtMinThresh'] = DirtData[1][1]
                    
                # Instantiate the Area, Num, and Sizes dictionaries
                ret[f]['Area'] = {}
                ret[f]['Num'] = {}
                ret[f]['Sizes'] = {}
                
                # Populate the Area dictionary with dirt area comparison data
                ret[f]['Area']['testAreaDirt'] = testArea
                ret[f]['Area']['stdAreaDirt'] = stdArea
                ret[f]['Area']['FractCorrect'] = round(AreafractCorrect,4)
                ret[f]['Area']['FractError'] = round(Areaerror,4)
                
                # Populate the Num dictionary with dirt number comparison data
                ret[f]['Num']['testNumDirt'] = testNum
                ret[f]['Num']['stdNumDirt'] = stdNum
                ret[f]['Num']['FractCorrect'] = round(NumfractCorrect,4)
                ret[f]['Num']['FractError'] = round(Numerror,4)
                
                # Populate the Sizes dictionary with dirt sizes comparison data
                try:
                    ret[f]['Sizes']['stdMaxSize'] = stdSizes.max()
                    ret[f]['Sizes']['testMaxSize'] = testSizes.max()
                    ret[f]['Sizes']['stdMinSize'] = stdSizes.min()
                    ret[f]['Sizes']['testMinSize'] = testSizes.min()
                    ret[f]['Sizes']['stdMeanSize'] = round(stdSizes.mean(),4)
                    ret[f]['Sizes']['testMeanSize'] = round(testSizes.mean(),4)
                    ret[f]['Sizes']['stdMedianSize'] = np.median(stdSizes)
                    ret[f]['Sizes']['testMedianSize'] = np.median(testSizes)
                except ValueError:
                    print f
                    print 
        else:
            print "Standard " +f+ "does not have a thresholds.txt file."
            print "Or " + pathtodir+" is not a directory."
    else:
        print "Finished"
    return ret
               
                

                
                
                