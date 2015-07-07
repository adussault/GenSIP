"""

Contains tools for performing operations on folders and many items in a folder.

"""

import os
import GenSIP.functions as fun
import cv2
from GenSIP.sandbox.Histograms import findMoPeakByImg

###################################################################################

###################################################################################

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

###################################################################################

###################################################################################

def DoForFolder(folder, function, *args,**kwargs):
    """
    Applies a function to all the items in a folder and returns a dictionary whose
    keys are the name of each file in the folder and whose values are the 
    results of applying the function to file. 
        Arguments:
            - folder: path to folder that contains all the files to be processed
            - function: function to apply to all the files in the folder
                NOTE: The first argument of this function must be a file path.
            - *args: arguments to pass on to the function in addition to the 
                     file path.
            - **kwargs: Keyword arguments to be passed on to the function as well
                     as keyword arguments for DoForFolder.
        Keyword Arguments:
            - filter: A function to apply to the folder contents to filter out any
                     undesirable items. 
    """
    
    contents = os.listdir(folder)
    filterfunction = kwargs.pop("filter",None)
    
    if filterfunction != None:
        contents = filterfunction(contents)
        
    ret = {}
    for sub in contents:
        path = folder + '/' + sub
        result = function(path,*args,**kwargs)
        name = os.path.splitext(sub)[0]
        ret[name] = result
    return ret
        
def DoImgForFolder(folderpath, function, *args,**kwargs):
    imgType = kwargs.get('imgType','tiff')
    subImgs = os.listdir(folderpath)
    subImgs[:] = [img for img in subImgs if img.startswith("sub_")]
    subImgs[:] = [img for img in subImgs if img.endswith(imgType)]
    subImgs[:] = [img for img in subImgs if not(img.endswith("000"+imgType))]
    subImgs[:] = [img for img in subImgs if not(img.startswith("sub_000"))]
    ret = []
    sbs = []
    for sub in subImgs:
        img = fun.loadImg(folderpath+'/'+sub)
        peak = function(img)
        sbs.append(sub)
        ret.append(peak)
    return ret,sbs
    
def findMoPeakByImgFolder(folderpath, imgType = ".tif", color=cv2.CV_LOAD_IMAGE_GRAYSCALE):
    subImgs = os.listdir(folderpath)
    subImgs[:] = [img for img in subImgs if img.startswith("sub_")]
    subImgs[:] = [img for img in subImgs if img.endswith(imgType)]
    subImgs[:] = [img for img in subImgs if not(img.endswith("000"+imgType))]
    subImgs[:] = [img for img in subImgs if not(img.startswith("sub_000"))]
    ret = []
    sbs = []
    for sub in subImgs:
        img = fun.loadImg(folderpath+'/'+sub)
        peak = findMoPeakByImg(img)
        sbs.append(sub)
        ret.append(peak)
    return ret,sbs