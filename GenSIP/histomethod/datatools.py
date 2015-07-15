""" 

Contains all the functions used by Histograms to analyze data on the plots. 
Unlike the histogram_tools module, this module has no functions specific only to
GenSIP data formats. 

"""
import numpy as np
from math import copysign

###################################################################################

###################################################################################
    
def easyThresh(img, MIN,MAX,Binary=True):
    """
    Does exactly the same thing as the threshold function, it just does it in 
    fewer lines of code, and doesn't have a Boolean option right now. 
    """
    ret = (MIN<=img)&(img<=MAX)
    if Binary:
        return ret
    else:
        return img[ret]
        
###################################################################################

###################################################################################
        
def atFWHM(histo, peakVal):
    """
    Returns the color indices at the full-width-half-maximum of a given peak
    on a histogram. Does not subtract the backround.
    """
    MAX = histo[peakVal]
    HalfMax = MAX/2
    colors = np.arange(0,histo.size,1)
    for i in colors[::-1][histo.size-peakVal:]:
        if histo[i]<HalfMax:
            HM1 = i
            break
    else:
        HM1=None
        #raise Exception("No first half-maximum.")
    
    for i in colors[peakVal:]:
        if histo[i]<HalfMax:
            HM2 = i
            break
    else:
        HM2=None
        #raise Exception("No second half-maximum.")
    return (HM1,HM2)

###################################################################################

###################################################################################

def nsmooth(a, i=4, n=3):
    """
    Runs smoothed over an array i times with a range n over which the running 
    mean is taken. 
    """
    ret = a.copy()
    for I in range(i):
        ret = smoothed(ret,n)
    return ret
         
###################################################################################

###################################################################################
 
def smoothed (a, n=3):
    """
    Wraps moving average in such a way that it returns an array of the same size
    as the input array. 
        Arguments:
         - a - a numpy nd array
        Key-word Arguments:
         - n - the number of points over which the average is calculated
    """
    mA = movingAverage(a)
    ret = np.zeros(a.shape)
    ret[0]=a[0]
    ret[255]=a[255]
    ret[1:255]=mA[:]
    return ret

###################################################################################

###################################################################################
 
def movingAverage(a,n=3):
    """
    Calculates the moving average at each point for an array of data and returns
    an array of the moving Averages. The length of that array is the length of the
    original array minus 2. 
    Arguments:
         - a - a numpy nd array
        Key-word Arguments:
         - n - the number of points over which the floating average is calculated
    
    """
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:]-ret[:-n]
    return ret[n-1:]/n

##################################################################################

###################################################################################
    
def runningMeanFast(x, N):
    """ 
    Another running mean function
        Arguments:
         - x - a numpy nd array
         - N - the number of points over which the floating average is calculated
    """
    return np.convolve(x, np.ones((N,))/N)[(N-1):]

###################################################################################

###################################################################################

def getMaxima(a, smoonum=4):
    """
    Returns all of the local maxima in a plot by first smoothing it 4 times 
    with a running mean and then looking for local maxima.
    """

    smoo = a.copy()
    for i in range(smoonum):
        smoo = smoothed(smoo)
    grad = np.gradient(smoo)
    Zx,Zy = getZeroes(grad)
    concav = np.gradient(grad)
    retX,retY = [],[]
    for x in Zx:
        if concav[x]<0:
            retX.append(x)
            retY.append(a[x])
    return retX,retY

###################################################################################

###################################################################################

def getMinima(a, smoonum=4):
    """
    Returns all of the local minima in a plot by first smoothing it 4 times 
    with a running mean and then looking for local minima.
    """
    smoo = a.copy()
    for i in range(smoonum):
        smoo = smoothed(smoo)
    grad = np.gradient(smoo)
    Zx,Zy = getZeroes(grad)
    concav = np.gradient(grad)
    retX,retY = [],[]
    for x in Zx:
        if concav[x]>0:
            retX.append(x)
            retY.append(a[x])
    return retX,retY

###################################################################################

###################################################################################
    
def getInflectionPoints(a, smoonum=4,sign='negative'):
    """
    Returns all of the negative or positive inflection points in a plot by first
    smoothing it 4 times with a running mean and then looking for points where 
    the second gradient is zero and the first gradient is either postive or negative.
    Inputs:
        a - an array of numbers
    Key-Word Arguments (kwargs):
        smoonum = 4 - number of times to smooth the array to remove noise. Default
            set to 4. 
        sign = 'negative' - specifies whether to look for positive or negative 
            inflection points. 
    """
    smoo = a.copy()
    for i in range(smoonum):
        smoo = smoothed(smoo)
    concav = 10*np.gradient(np.gradient(smoo))
    Zx,Zy = getZeroes(concav)
    aberr = np.gradient(concav)
    retX,retY = [],[]

    for x in Zx:
        if sign=='negative':
            if aberr[x]<0:
                retX.append(x)
                retY.append(a[x])
        elif sign=='positive':
             if aberr[x]>0:
                retX.append(x)
                retY.append(a[x])
        else:
            retX.append(x)
            retY.append(a[x])
    return retX,retY

###################################################################################

###################################################################################

def getZeroes(a):
    """Returns the index of the points at which an array crosses zero given the array"""
    x,y = [],[]
    for i in range(a[0:-1].size):
        if a[i] == 0:
            x.append(i)
            y.append(a[i])
        elif copysign(a[i],a[i+1])!=a[i]:
            x.append(i)
            y.append(a[i])
    return x,y
                     