"""
This module contains the functions used to deal with a large panorama image such 
as those used in analyzing the big foil scans. 

"""

import numpy as np
import os
import GenSIP_v4.functions as fun
import cv2
import GenSIP_v4.kuwahara as K
from scipy import misc

###################################################################################

###################################################################################

def splitImage(image, numParts, path = "InputPicts/FoilScans/Q1",name=""):
    """
    Takes an image and divides it up into a number of sub-images specified by 
    numParts. writes the output to a folder in the path folder and names it 
    "sub_imgs_"+name.
        Inputs:
         - image - the image to be divided
         - numParts - the number of sub-Images to be produced. Must be a perfect
            square.
        Key-word Arguments:
         - path - path to the folder to contain the subimages folder
         - name - name of the images produced
    """
    
    perSide = int(np.sqrt(numParts))
    height = image.shape[0]
    width = image.shape[1]
    
    h_unit = int(height/perSide)
    w_unit = int(width/perSide)
    h_rem = int(height%perSide)
    w_rem = int(width%perSide)
    start_col = 0
    start_row = 0
    
    # Create output folder
    outPath = path+"/sub_imgs_"+name
    if not(os.path.exists(outPath)):
        os.makedirs(outPath)
    
    
    for r in range(perSide): # Column
        #print "Starting Row " + str(r)
        stop_row = h_unit + r*h_unit
        
        if r == perSide-1:
            stop_row += h_rem - 1
        
        for c in range(perSide): # Row
            #print "Row: " + str(r) +" Column: " + str(c)
            stop_col = w_unit + c*w_unit
            if c == 0:
                start_col = 0
            if c == perSide-1:
                stop_col += w_rem - 1
            subImage = image[start_row:stop_row,start_col:stop_col]
            '''
            print "Subimage Made. Writing Image..."
            print "Start_row: " + str(start_row)
            print "Stop_row: " + str(stop_row)
            print "Start_col: " + str(start_col)
            print "Stop_col: " + str(stop_col)
            print "Image shape: " + str(subImage.shape)
            '''
            cv2.imwrite(str(outPath+"/sub_"+str(r).zfill(3) +"_"+str(c).zfill(3)+".tif"),subImage)
        
            start_col = stop_col+1
            
        start_row = stop_row+1

###################################################################################

###################################################################################

def stitchImage(folderpath, imgType = ".png", color=cv2.CV_LOAD_IMAGE_GRAYSCALE):
    """
    This method takes a folder of sub-images and stitches them back together and
    writes the resulting image to a large montage image named "montage.png" in 
    the given folder.
    
    Input Arguments:
        - folderpath - path to the folder containing the sub-images
    Key-Word Arguments:
        - imgType - type of image file contained in subfolder (i.e. .jpeg, .png,
          .tif, etc.) it is important to include the '.' at the beginning. Default
          set to ".png".
        - color - specifies whether the images are to be loaded in color or 
          grayscale. 

    """

    subImgs = os.listdir(folderpath)
    subImgs.sort()
    montageHeight = 0
    montageWidth = 0
    
    # Make sure that the folder of sub-images only contains the sub images to be
    # stitched together.
    """
    for img in subImgs:
        if not(img.startswith("sub_")):
            subImgs.remove(img)
            print img
        elif not(img.endswith(".png")):
            subImgs.remove(img)
    """
    
    subImgs[:] = [img for img in subImgs if img.startswith("sub_")]
    subImgs[:] = [img for img in subImgs if img.endswith(imgType)]
    
    if len(subImgs) == 0:
        raise Exception("You may have specified the wrong image Type.")

    def getIndex(name):
        index = name.strip("sub_")
        index = index.strip(imgType)
        index = index.rsplit("_")
        index = [int(index[0]),int(index[1])]
        return index
        
    for img in subImgs:
        index = getIndex(img)
        sub = fun.loadImg(folderpath+"/"+img,0)
        if index[1]==0:
            montageHeight += sub.shape[0]
        if index[0]==0:
            montageWidth += sub.shape[1]
        del(sub)
    #return (montageHeight,montageWidth)
    # See if image is meant to be loaded in color or not. If so, make 
    # the montage a color image (a 3D array with the 3rd dimension as RGB)
    if len(fun.loadImg(folderpath+"/"+subImgs[0],color).shape)==3:
        is3D = True
        montage = np.zeros((montageHeight,montageWidth,3))
    else:
        is3D = False
        montage = np.zeros((montageHeight,montageWidth))
        
    # Initialize iterating variables for the for loop:
    lastindex = [0,0]
    lastSubW = 0
    lastSubH = 0
    height = 0
    width = 0
    
    for img in subImgs:
        # index gives the tuple index of the sub-image in the montage matrix
        index = getIndex(img)
        
        #list index gives the index of the sub-image in the list of subImgs
        """
        listIndex = subImgs.index(img)
        
        if listIndex<len(subImgs)-1:
            nextIndex = getIndex(subImgs[listIndex+1])
        else:
            nextIndex = (0,0)
        """
        
        # Load the sub image
        sub = fun.loadImg(folderpath+"/"+img,color)
        # Assign the height and width of the sub-image to subH and subW
        # Check to see if the image is in color when loaded. If so, convert
        # from BGR to RGB since Opencv by default loads the image as BGR.
        if is3D:
            subH,subW,rgb = sub.shape
            sub=cv2.cvtColor(sub,cv2.COLOR_BGR2RGB)
        else:
            subH,subW = sub.shape
        
        # Set the variables for the indices in the montage over which the 
        # sub-image will be inserted. These variables are set according to 
        # the relationship between the current sub-image index and the last
        # sub-image index, to see if the current sub-image starts a new row.
        if index[0] == 0:
        # if the sub-image is part of the first row
            lastSubH = 0
            height = subH
        elif index[0] == lastindex[0]+1:
        # if the sub-image is the start of a new row but is not
        # in the first row. 
            lastSubH = height
            height += subH
            
        if index[1] == 0:
        # if the sub-image is the start of a new row
            lastSubW = 0
            width = subW
            
        elif index[1] == lastindex[1]+1:
        # if the sub-image continues a row
            lastSubW = width
            width += subW
        print img
        
        # Insert the sub-image into the montage:
        montage[lastSubH:height,lastSubW:width]=sub[:,:]
        # Iterate the variable lastindex:
        lastindex = index
        
    # Write the montage image to the folder containing
    cv2.imwrite(folderpath+"/montage.png",montage,[cv2.cv.CV_IMWRITE_PNG_COMPRESSION,6])
    
###################################################################################

###################################################################################

def makeManyPosters(foldername,foilname="40360_2", Quarter="Q1",Mask=0, kern=6,\
KuSize=9,Gaus1=3,Gaus2=11,rsize=.1,Kuw_only=False, ExcludeDirt=True):
    """
    Runs BIGmakeposter on a set of images in a folder and saves them to a new 
    poster folder.
    """
    images = os.listdir(foldername)
    
    if Kuw_only:
        maps="/KuwaharaMaps"
    else:
        maps = "/PosterMaps"
    
    outFolder = 'Output/Output_'+foilname+"_"+Quarter
    if not os.path.exists(outFolder):
        os.makedirs(outFolder)
    if not os.path.exists(outFolder+maps):
        os.makedirs(outFolder+maps)
    
        
    for image in images:
        img = fun.loadImg(foldername+'/'+image)
        root,ext = os.path.splitext(image)
        poster = bigPostPreProc(img,Mask,kern,KuSize,Gaus1,Gaus2,rsize,Kuw_only,ExcludeDirt)
        cv2.imwrite(outFolder+maps+'/'+root+".png",poster, \
        [cv2.cv.CV_IMWRITE_PNG_COMPRESSION,6])
        print image
  
###################################################################################

###################################################################################
     
def bigPostPreProc(image,**kwargs):
    """
    This method takes the image of the foil and creates a smoothed Kuwahara image
    used to make the poster for regional thresholding.
        Inputs:
         - image - the image to perform preprocessing on.
        Key-Word Arguments:
             - Mask = 0 - Assign the Mask here
             - KuSize = 17 - Size of Kuwahara Filter
             - Gaus1 = 3 - Size of first Gaussian Blur
             - Gaus2 = 11 - Size of second Gaussian Blur
             - rsize = .1 - Resize value
             - Kuw_only = False - Option to only return the Kuwahara filtered image
             - ExcludeDirt = True - Option to Exclude dirt 

    """
    # The keyword 
    Mask = kwargs.get("Mask",0) # Assign the Mask here
    KuSize = kwargs.get("KuSize",17) # Size of Kuwahara Filter
    Gaus1 = kwargs.get("Gaus1",3) # Size of first Gaussian Blur
    Gaus2 = kwargs.get("Gaus2",11) # Size of second Gaussian Blur
    rsize = kwargs.get("rsize",.1) # Resize value
    Kuw_only = kwargs.get("Kuw_only",False) # Option to only return the Kuwahara filtered image
    ExcludeDirt = kwargs.get("ExcludeDirt",True) # Option to Exclude dirt 

    img = np.copy(image)
    averageColor = int(np.average(img))
    if ExcludeDirt:
        img[img<=40]=averageColor
    if type(Mask)==np.ndarray and Mask.shape==image.shape:
        invMsk = np.bitwise_not(Mask)
        invMsk = int(np.average(image))*(invMsk/255)
        img = cv2.add(image,invMsk)
        img[img>np.max(image)] = int(np.average(image))
        
    rsz = misc.imresize(img,rsize,interp='bicubic')
    if Kuw_only:
        return rsz
    gr = cv2.GaussianBlur(rsz, (Gaus1,Gaus1),0)
    kgr = K.Kuwahara(gr,KuSize)
    kgr = kgr.astype(np.uint8) # Make sure the Kuwahara image is uint8 so it doesn't scale
    rkgr = misc.imresize(kgr,(image.shape),interp='bicubic')
    grkgr = cv2.GaussianBlur(rkgr, (Gaus2,Gaus2),0)
    if Kuw_only:
        return rkgr
    else:
        return grkgr

###################################################################################

###################################################################################
       
def bigPosterfy(image,k_size=6):
    """
    Takes a gray image and sets all values within a given range to a single value
    this way we can divide up regions into high Exposure areas, primarily Pt, 
    primarily Mo, dirt, or black. The only difference between bigPosterfy and 
    posterfy found in GenSIP_v4.functions is the additional 'highEx' region. This
    is added because the large foils have more fluctuation in topography, and 
    hence various levels of contrast and exposure. 
    
    Kwargs:
        - k_size - size of the kernal used in the final morphological opening
            step. 
    """
    image_copy=image.copy()
    image_copy=image_copy.astype(np.uint8)
    blk=(image_copy<=4)
    pleat=(5<=image_copy)&(image_copy<=40)
    darkMo=(41<=image_copy)&(image_copy<=139)
    Mo=(140<=image_copy)&(image_copy<=195)
    highEx=(196<=image_copy)&(image_copy<=215)
    Plat=(216<=image_copy)
    image_copy[blk]=0
    image_copy[pleat]=50
    image_copy[darkMo]=85
    image_copy[Mo]=150
    image_copy[highEx]=200
    image_copy[Plat]=255
    #image_copy = maskEdge(image_copy,thickness=60)
    #Do a morphological opening step to eliminate the rough edges
    # I prefer opening over closing because it is more important to catch
    # all of the dark areas, as the Pt areas will have pretty sraightforward
    # threshold results
    kernel = np.ones((k_size,k_size))
    #image_copy = cv2.morphologyEx(image_copy, cv2.MORPH_OPEN, kernel)
    return image_copy