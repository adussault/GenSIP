"""
Contains thresholding, masking, and flooding operations. Called 'binaryops' 
because these functions all produce or deal with binary images.

Functions include:
   floodByThresh
   floodBySeed 
   FloodStep
   threshold
   easythreshold
"""

import numpy as np

def floodByThresh(img, min0, max0, d, upBound, lowBound, growthMin=0, growthMax=10000):
    """ 
    A Floodfill algorithm that takes an image and, starting with all the pixels 
    in a certain range, expands those regions to take in pixels with a similar 
    graylevel. It is meant for filling out dirt paricles when a portion of the 
    dirt particle is known to be dirt and another portion is not certainly dirt.
        Inputs:
         - img - the image to be flooded
         - min0 - the minimum value for the initial threshold of the image
         - max0 - the maximum value for the initial threshold of the image
         - d - the maximum difference in pixel value that will result in a pixel
            being 'flooded' in a step.
         - upBound - the uppermost graylevel that will be considered able to be 
            flood-filled. Should be greater than max0.
         - lowBound - the lowermost graylevel that will be considered able to be 
            flood-filled. Should be less than min0.
         
        Key-word Arguments:
         - growthMin - the minimum number of pixels that the Flooder is allowed 
            to grow in one step. Once the growth in a step is <= growthMin, the 
            Flooder stops. 
         - growthMax - the maximum number of pixels that the Flooder is allowed 
            to grow in one step. Once the growth in a step is >= growthMax, the 
            Flooder stops. Prevents the flooder from incoporating all 

    """
    if growthMin < 0:
        raise Exception("Minimum growth level must be 0 or positive.")
    ret = threshold(img, min0,max0,Binary=True,Bool=True)
    growth = growthMin+1
    i = 1
    while (growthMin<growth<growthMax):
        print "Step #" + str(i) + ":"
        ret, growth = FloodStep(img, ret, d, upBound, lowBound, growthMin, growthMax)
        i += 1
        if i == 100:
            break
    return ret

def floodBySeed(img, seed, d, upBound, lowBound, growthMin=0, growthMax=10000, verbose=False):
    """ 
    A Floodfill algorithm that takes an image and, starting with all the pixels 
    in a certain range, expands those regions to take in pixels with a similar 
    graylevel. It is meant for filling out dirt paricles when a portion of the 
    dirt particle is known to be dirt and another portion is not certainly dirt.
        Inputs:
         - img - the image to be flooded
         - seed - a binary image that sets the starting point for the flooder
         - d - the maximum difference in pixel value that will result in a pixel
            being 'flooded' in a step.
         - upBound - the uppermost graylevel that will be considered able to be 
            flood-filled. Should be greater than max0.
         - lowBound - the lowermost graylevel that will be considered able to be 
            flood-filled. Should be less than min0.
         
        Key-word Arguments:
         - growthMin - the minimum number of pixels that the Flooder is allowed 
            to grow in one step. Once the growth in a step is <= growthMin, the 
            Flooder stops. 
         - growthMax - the maximum number of pixels that the Flooder is allowed 
            to grow in one step. Once the growth in a step is >= growthMax, the 
            Flooder stops. Prevents the flooder from incoporating all 

    """
    if growthMin < 0:
        raise Exception("Minimum growth level must be 0 or positive.")
    ret = seed.copy()
    growth = growthMin+1
    i = 1
    while (growthMin<growth<growthMax):
        if verbose:
            print "Step #" + str(i) + ":"
        ret, growth = FloodStep(img, ret, d, upBound, lowBound, \
        growthMin=growthMin, growthMax=growthMax, verbose=verbose)
        i += 1
        if i == 1000:
            break
    return ret

def FloodStep(img, seed, d, upBound, lowBound, growthMin=0, growthMax=1000, verbose=False):
    """
    Does one step of the flooding process by looking at the nearest neighbors of
    the seed pixels (pixels that are already flooded) and adding those that have 
    pixel values within d to the flooded pixels.
        Inputs:
         - img - image over which the flood step is performed
         - seed - a binary (black & white) image where all the white pixels are
            the initial, already flooded, 
    """
    seed = seed.astype(np.bool_)
    ret = seed.copy()
    imask,jmask = np.where(seed)#old: ((img>=min0)&(img<=max0))
    imask,jmask = imask.copy(),jmask.copy() #Make writable
    #ret = np.zeros(img.shape).astype(np.bool_)
    #ret[imask,jmask] = True
    ni = (-1,0,1)
    nj = (-1,0,1)
    growth = 0
    for di in ni:
        for dj in nj:
            if not((di==dj) and (di==0)):
                iself = imask
                jself = jmask
                ineigh = imask + di
                jneigh = jmask + dj
                # Check the Borders:
                ineigh,jneigh,iself,jself = \
                ineigh[(ineigh>=0)&(jneigh>=0)&(ineigh<img.shape[0])&(jneigh<img.shape[1])],\
                jneigh[(ineigh>=0)&(jneigh>=0)&(ineigh<img.shape[0])&(jneigh<img.shape[1])],\
                iself[(ineigh>=0)&(jneigh>=0)&(ineigh<img.shape[0])&(jneigh<img.shape[1])],\
                jself[(ineigh>=0)&(jneigh>=0)&(ineigh<img.shape[0])&(jneigh<img.shape[1])]
                
                # Eliminate pixels that are already counted
                ineigh,jneigh,iself,jself = \
                ineigh[seed[ineigh,jneigh]==False],\
                jneigh[seed[ineigh,jneigh]==False],\
                iself[seed[ineigh,jneigh]==False],\
                jself[seed[ineigh,jneigh]==False]
                #print str(img[ineigh,jneigh])
                
                # Compare values of pixels and neighbor pixels
                selfVals = img[iself,jself]
                neighVals = img[ineigh,jneigh]
                
                # Calculate the absolute difference between the pixel and its neighbors.
                # First convert to integer to prevent overflow of the 8-bit 
                # unsigned pixel values, then convert back to uint8
                absdiff = np.absolute(selfVals.astype(np.int)-neighVals.astype(np.int)).astype(np.uint8)

                # Get i,j indices of qualifying neighbors
                iquals = ineigh[(absdiff<=d)&(neighVals<=upBound)&(neighVals>=lowBound)]
                jquals = jneigh[(absdiff<=d)&(neighVals<=upBound)&(neighVals>=lowBound)]
                # increase growth by the number of qualifying neighbors
                growth += iquals.size
                
                #print str(img[iquals,jquals])
                
                # set qualifying neighbors in return image to True
                ret[iquals,jquals] = True
    if verbose:           
        print "  Growth: " + str(growth)
        print "  Less than max? " + str(growth<growthMax)
        print "  More than min? " + str(growth>growthMin)
    return ret, growth
    
    '''
    imask,jmask = np.where((img>=min0)&(img<=max0))
    imask,jmask = imask.copy(),jmask.copy()
    if Binary:
        ret = np.zeros(img.shape).astype(np.bool_)
        ret[imask,jmask] = True
    else:
        ret = np.zeros(img.shape).astype(img.dtype)
        ret[imask,jmask] = img[imask,jmask]
    return ret
    '''
def threshold(img, MIN,MAX=255,Binary=True,Bool=False):
    '''
    A home made simple thresholding function. Not as robust as the one in opencv
    since it only has one straightforward thresholding method, and doesn't have
    all the options as in opencv. But it works.
    Returns a thresholded image of either type np.uint8 or bo
    Inputs are:
        Arguments:
            
        Key-word Arguments:
         - MAX - The maximum value of the threshold. Default set to 255.
         - Binary - Option to output a black and white uint8 image, with black 
            set to 0 and white set to 255
        - Boolean - Option to set the output image to np.bool_ type, which means
            that black will be 0 (or False) and white would be 1 (or True).
    
    A nice student project would be to get this thresholding method working for 
    RGB images. It is low priority though since we deal with grayscale images.
    
    '''
    ret = img.copy()
    if len(ret.shape)==2:
        ret[ret<MIN] = 0
        ret[ret>MAX] = 0
        if Binary:
            if Bool:
                ret=ret.astype(np.bool_)
                ret[ret!=0] = True
            else:
                ret[ret!=0] = 255
        return ret
    elif len(ret.shape)==3 and ret[2]==3:
        raise Exception("I only deal with grayscale images for now")