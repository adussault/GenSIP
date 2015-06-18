# -*- coding: utf-8 -*-
"""
This module contains functions used to analyze the precision and accuracy of the
GenSIP analysis tools, and set boundaries on the error in the results of the GenSIP
analysis. 
"""
from GenSIP.functions import *
from GenSIP.molybdenum import getFoilArea, calcPtArea
from GenSIP.dirt import amountDirt
from matplotlib import pyplot as plt
from socket import gethostname
from time import time
import pickle # Used to serialize the "Data" dictionary.
import numpy as np

def threshPrec(sss, Domain=120, step=2, MoDirt="Mo",RUN=("P","D","M","Pt"),Pstart=70,Dstart=100,Mstart=160,Ptstart=160):
    """
    This method is used to examine the error of analyzeMo and analyzeDirt
    due to variations in the threshold level. It produces a series of graphs 
    depicting the variation of the calculated molybdenum loss or dirt loss with 
    the threshold levels in the regionalThresh.
        Inputs:
            - sss: The Sample Set String of input pictures to examine, just as in
            - Domain: An integer that gives the length of the threshold domain 
                over which to plot the output values. (Sorry that is terribly 
                explained. Its the length of the x-axis in the output plot.) 
            - step: Step length over the Domain. Default set to 2
            - MoDirt: string variable designating whether to examine analyseMo or
                analyzeDirt. Allowable values are:
                    - Moly Analysis: "Mo","Moly","moly","molybdenum","Molybdenum","M","m"
                    - Dirt Analysis: "Dirt","dirt","D","d"
            - RUN: a tuple describing which regions you want to vary the threshold
                value over. Regions that are not in the RUN tuple will be kept 
                constant at their start value (below). 
            - Pstart, Dstart, Mstart, Ptstart: The values at which to start the 
                threshold domain. 
            
            
    For Mo analysis:
        Pstart=90,Dstart=120,Mstart=180,Ptstart=180
    For dirt analysis:
        Pstart=8,Dstart=28,Mstart=55,Ptstart=60
    
    """
    t1 = time()
    # Create arrays for the plot:
    numpoints = int(Domain/step)
    
    # The x-axis will start at Mstart, since that is typically the most dominant
    # region of the foil:
    x = np.linspace(Mstart,Mstart+Domain, numpoints).astype(np.uint16)
    
    # Set up domain in different regions:
    if "P" in RUN:
        P = np.linspace(Pstart,Pstart+Domain, numpoints).astype(np.uint16)
        P[P>255] = 255
    else:
        P = np.zeros((numpoints)).astype(np.uint16)
        P[:] = Pstart
        
    if "D" in RUN:
        D = np.linspace(Dstart,Dstart+Domain, numpoints).astype(np.uint16)
        D[D>255] = 255
    else:
        D = np.zeros((numpoints)).astype(np.uint16)
        D[:] = Dstart
        
    if "M" in RUN:
        M = np.linspace(Mstart,Mstart+Domain, numpoints).astype(np.uint16)
        M[M>255] = 255
    else:
        M = np.zeros((numpoints)).astype(np.uint16)
        M[:] = Mstart
        
    if "Pt" in RUN:
        Pt = np.linspace(Ptstart,Ptstart+Domain, numpoints).astype(np.uint16)
        Pt[Pt>255] = 255
    else:
        Pt = np.zeros((numpoints)).astype(np.uint16)
        Pt[:] = Ptstart
        
        
    # Generate plot subtitle:
    subtitle = "varies "
    Startstring = "from: "
    starts = {"P":Pstart,"D":Dstart,"M":Mstart,"Pt":Ptstart}
    for r in RUN:
        subtitle = r +", " + subtitle
        Startstring = Startstring+str(starts.pop(r))+", "
    End = "Constant: "
    for i in starts:
        End = End + i + " at "+str(starts[i])+", "
    Startstring = Startstring[:-2]+". " + End[:-2]+"."
    subtitle = subtitle + Startstring
        
        
    #print "All runners right size?  " + str(P.size==D.size==M.size==Pt.size)
    
    # Declare a list of acceptable file types: tif and jpg
    filetypes = ['.tif', '.jpg', '.jpeg','.tiff']
    # Make list of files in the corresponding Before and After folders
    befpics = sorted(os.listdir('InputPicts/Before/Before_'+sss))
    aftpics = sorted(os.listdir('InputPicts/After/After_'+sss))
    
    # Make output folder if does not exist
    if not os.path.exists('Output/Output_'+sss):
        os.makedirs('Output/Output_'+sss)
    if not os.path.exists('Output/Output_'+sss+'/Plots/'+'Threshold_'+MoDirt):
        os.makedirs('Output/Output_'+sss+'/Plots/'+'Threshold_'+MoDirt)
    
    allowableMo = ("Mo","Moly","moly","molybdenum","Molybdenum","M","m")
    allowableDirt = ("Dirt","dirt","D","d")
    
    # Create a Dictionary to store the data:
    datestring = getDateString()
    host = os.path.splitext(gethostname())[0]
    Version = getGenSIPVersion()
    
    Data = {
    "Date":datestring,"Host":host, 
    "Version":Version,"SSS":sss,"Analysis":MoDirt,"Domain":x,
    "SETTINGS":{"RUN":RUN,"Pstart":Pstart,"Dstart":Dstart,"Mstart":Mstart,"Ptstart":Ptstart},
    "DATA":{}
    }
    Data["Domain"] = x
    
    
    for fn in befpics:
        filename,exten = os.path.splitext(fn)
        befpath = 'InputPicts/Before/Before_'+sss+'/'+fn
        if exten in filetypes:
            # Get number of current foil
            foilnum = filename.split()[0]
            aftfoil = foilnum + ' after clean' + exten
            aftpath = 'InputPicts/After/After_'+sss+'/'+aftfoil
            # Check to see if the corresponding picture is in the after folder
            if aftfoil in aftpics:
                
                y1 = np.zeros((numpoints))
                y2 = np.zeros((numpoints))
                yBEF = np.zeros((numpoints))
                yAFT = np.zeros((numpoints))
                
                if MoDirt in allowableMo:
                    for i in range(numpoints):
                        ret = CalMoComp(befpath, aftpath, P=P[i],D=D[i],M=M[i],Pt=Pt[i],I=str(i))
                        Ptbef = ret[0]
                        Ptaft = ret[1]
                        #PerPtLoss = (Ptbef-Ptaft)/Ptbef
                        massLoss = ret[3]
                        PercLoss = ret[4]
                        y1[i] = massLoss
                        y2[i] = PercLoss
                        yBEF[i] = Ptbef
                        yAFT[i] = Ptaft
                    
                    # MAKE PLOT 
                    
                    plt.figure()
                    plt.suptitle(foilnum+" "+MoDirt+" Analysis")
                    plt.title(subtitle,fontsize=14)
                    plt.plot(x, y1, "-ro")
                    plt.plot(x, y2, "-bo")
                    plt.legend(["mass Mo Lost","% Mo Loss"],loc=0)
                    plt.xlabel("Threshold Values")
                    plt.ylabel("Moly Loss")
                    #Set Max and min of the y axis:
                    minY = min(np.min(y1),np.min(y2)) - 2
                    maxY = max(np.max(y1),np.max(y2)) + 5
                    if minY >= 0:
                        minY = 0
                    plt.axis([x[0]-2,x[numpoints-1]+2,minY, maxY])
                    plt.savefig('Output/Output_'+sss+'/Plots/'+'Threshold_'+MoDirt+"/"+foilnum+" Mo Threshold plot.png")
                    print "Moly plot for "+foilnum+" finished."
                    plt.close()
                    
                    plt.figure()
                    plt.suptitle(foilnum+" "+MoDirt+" Analysis")
                    plt.title(subtitle,fontsize=14)
                    plt.plot(x, yBEF, "-ro")
                    plt.plot(x, yAFT, "-bo")
                    plt.legend(["Pt count before","Pt count after"],loc=0)
                    plt.xlabel("Threshold Values")
                    plt.ylabel("Pt count")
                    #Set Max and min of the y axis:
                    minY = min(np.min(y1),np.min(yBEF)) - 2
                    maxY = max(np.max(y1),np.max(yAFT)) + 5
                    if minY >= 0:
                        minY = 0
                    plt.axis([x[0]-2,x[numpoints-1]+2,minY, maxY])
                    plt.savefig('Output/Output_'+sss+'/Plots/'+'Threshold_'+MoDirt+"/"+foilnum+" Pt area plot.png")
                    print "Pt area plot for "+foilnum+" finished."
                    plt.close()
                    
                    plt.figure()
                    plt.suptitle(foilnum+" "+MoDirt+" Derivative")
                    plt.title(subtitle,fontsize=14)
                    grad1 = np.gradient(y1, step)
                    grad2 = np.gradient(y2, step)
                    plt.plot(x, grad1, "-ro", label = "Deriv Dirt loss by count")
                    plt.plot(x, grad2, "-bo", label="Deriv Dirt loss by area")
                    plt.legend(["Deriv mass Mo Lost","Deriv % Mo Loss"],loc=0)
                    plt.xlabel("Threshold Values")
                    plt.ylabel("dy/dx")
                    #Set Max and min of the y axis:
                    minY = min(np.min(grad1),np.min(grad2)) - 2
                    maxY = max(np.max(grad1),np.max(grad2)) + 5
                    if minY >= 0:
                        minY = 0
                    plt.axis([x[0]-2,x[numpoints-1]+2,minY, maxY])
                    plt.savefig('Output/Output_'+sss+'/Plots/'+'Threshold_'+MoDirt+"/"+foilnum+" Mo Gradient plot.png")
                    plt.close()
                    
                    Data["DATA"][foilnum] = {"massMoLoss":y1,"% MoLoss":y2,"massGrad":grad1,"percGrad":grad2}
                    
                    print "Moly gradient plot for "+foilnum+" finished."
                    
                elif MoDirt in allowableDirt:
                    for i in range(numpoints):
                        ret = CaldirtComp(befpath, aftpath,P=P[i],D=D[i],M=M[i],Pt=Pt[i])
                        numBef = ret[0]
                        numAft = ret[1]
                        areaBef = ret[2]
                        areaAft = ret[3]
                        percNum = 100*(numBef-numAft)/numBef
                        percArea = 100*(areaBef-areaAft)/areaBef
                        y1[i] = percNum
                        y2[i] = percArea
                    print "Dirt plot for "+foilnum+" finished."
                    
                    plt.figure()
                    plt.suptitle(foilnum+" "+MoDirt+" Analysis")
                    plt.title(subtitle,fontsize=14)
                    plt.plot(x, y1, "-ro", label = "Dirt loss by count")
                    plt.plot(x, y2, "-bo", label="Dirt loss by area")
                    plt.legend(["Dirt loss by count","Dirt loss by area"],loc=0)
                    plt.xlabel("Threshold Values")
                    plt.ylabel("Percent")
                    #Set Max and min of the y axis:
                    minY = min(np.min(y1),np.min(y2)) - 2
                    maxY = max(np.max(y1),np.max(y2)) + 5
                    if minY >= 0:
                        minY = 0
                    plt.axis([x[0]-2,x[numpoints-1]+2,minY, maxY])
                    plt.savefig('Output/Output_'+sss+'/Plots/'+'Threshold_'+MoDirt+"/"+foilnum+" Dirt Threshold plot.png")
                    plt.close()
                    
                    plt.figure()
                    plt.suptitle(foilnum+" "+MoDirt+" Derivative")
                    plt.title(subtitle,fontsize=14)
                    grad1 = np.gradient(y1, step)
                    grad2 = np.gradient(y2, step)
                    plt.plot(x, grad1, "-ro")
                    plt.plot(x, grad2, "-bo")
                    plt.legend(["Deriv Dirt loss by count","Deriv Dirt loss by area"],loc=0)
                    plt.xlabel("Threshold Values")
                    plt.ylabel("dy/dx")
                    #Set Max and min of the y axis:
                    minY = min(np.min(grad1),np.min(grad2)) - 2
                    maxY = max(np.max(grad1),np.max(grad2)) + 5
                    if minY >= 0:
                        minY = 0
                    plt.axis([x[0]-2,x[numpoints-1]+2,minY, maxY])
                    plt.savefig('Output/Output_'+sss+'/Plots/'+'Threshold_'+MoDirt+"/"+foilnum+" Dirt Gradient plot.png")
                    plt.close()
                    
                    plt.figure()
                    plt.suptitle(foilnum+" "+MoDirt+" Analysis")
                    plt.title(subtitle,fontsize=14)
                    plt.plot(x, y1, "-ro", label = "Dirt loss by count")
                    plt.plot(x, y2, "-bo", label="Dirt loss by area")
                    plt.legend(["Dirt loss by count","Dirt loss by area"],loc=0)
                    plt.xlabel("Threshold Values")
                    plt.ylabel("Percent")
                    #Set Max and min of the y axis:
                    minY = min(np.min(y1),np.min(y2)) - 2
                    maxY = max(np.max(y1),np.max(y2)) + 5
                    if minY >= 0:
                        minY = 0
                    plt.axis([x[0]-2,x[numpoints-1]+2,minY, maxY])
                    plt.savefig('Output/Output_'+sss+'/Plots/'+'Threshold_'+MoDirt+"/"+foilnum+" Dirt Threshold plot.png")
                    plt.close()
                    
                else:
                    raise Exception("\n MoDirt must be one of the following: \n"+\
                    str(allowableMo)+"\n - - - - - or - - - - - \n"+str(allowableDirt))
            else: print "No after image for foil "+foilnum
        else:
            print "Not a picture."
    # Get the run time of the function and add it to the dictionary "Data":
    t2 = time()-t1
    Data["TimeEllapse"] = convTime(t2)
    # Serialize the "Data" dictionary and store it in a pickle file with the other images. 
    with open('Output/Output_'+sss+'/Plots/'+'Threshold_Data_'+MoDirt+str(Domain),'wb') as f:
        pickle.dump(Data,f)
    
    
        
        
def compFoil():
    for fn in befpics:
        filename,exten = os.path.splitext(fn)
        befpath = 'InputPicts/Before/Before_'+sss+'/'+fn
        if exten in filetypes:
            # Get number of current foil
            foilnum = filename.split()[0]
            aftfoil = foilnum + ' after clean' + exten
            aftpath = 'InputPicts/After/After_'+sss+'/'+aftfoil
            # Check to see if the corresponding picture is in the after folder
            if aftfoil in aftpics:
        # Compare the two images:
                dirtVals, dirtPicts = dirtComp(befpath,aftpath)   
                (numbf,numaf,areabf,areaaf) = dirtVals
                (threshedbf, threshedaf) = dirtPicts
                # Calculate difference in area
                perDirtLoss = 100*(areaaf-areabf)/areabf
        	
        	# Write the output to a new line in the csv file
        	dirtwriter.writerow([foilnum, numbf, numaf, areabf, areaaf, perDirtLoss])
                # Create the output image
                # outfoil = compareMap(picts)
                # Save outfoil to Output folder
                # cv2.imwrite('Output/Output_'+sss+'/'+foilnum+' compare'+exten, outfoil)
                # Save the threshed images for analysis
                cv2.imwrite('Output/Output_'+sss+'/DirtMaps/'+foilnum+' before_threshed.png',\
                threshedbf, [cv2.cv.CV_IMWRITE_PNG_COMPRESSION,6])
                cv2.imwrite('Output/Output_'+sss+'/DirtMaps/'+foilnum+' after_threshed.png',\
                threshedaf, [cv2.cv.CV_IMWRITE_PNG_COMPRESSION,6])
                
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
	   print filename+" is not a picture."
        
        
        
def CalMoComp (before, after, P=90,D=120,M=180,Pt=180,res=1,I=""):
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
    bf = loadImg(before, cv2.CV_LOAD_IMAGE_GRAYSCALE)
    af = loadImg(after, cv2.CV_LOAD_IMAGE_GRAYSCALE)

    # Generate binary thresholds for Platinum:
    Ptbef = CALisolatePt(bf,P=P,D=D,M=M,Pt=Pt)
    Ptaft = CALisolatePt(af,P=P,D=D,M=M,Pt=Pt)

    # Approximate the percent of molybdenum lost:
    # Will take the average 
    foilareabef = getFoilArea(bf,res)
    foilareaaft = getFoilArea(af,res)

    # Calculate area of exposed platinum, adjust for resolution:
    PtAreaBef = calcPtArea(Ptbef)*res
    PtAreaAft = calcPtArea(Ptaft)*res
	
    #print "Total foil area before: %s" %str(foilareabef)
    #print "Total foil area after:  %s" %str(foilareaaft)

    # Calculate difference in area of exposed Pt:
    # Normalize values to the area of the before image
    befaftRatio = foilareabef/foilareaaft
    areaLoss = PtAreaAft*befaftRatio-PtAreaBef
    MolyBef = foilareabef-PtAreaBef
    MolyAft = foilareaaft-PtAreaAft
    #print "Mo area before: %s" %str(MolyBef)
    #print "Mo area after:  %s" % str(MolyAft)
    pcntMolyLoss = 1-(MolyAft/MolyBef)*(foilareabef/foilareaaft)
    pcntMolyLoss = round(pcntMolyLoss*100,1)
    # Now approximate molybdenum loss:
    # assuming 1pixel == 1 micron
    # Moly thickness ~300nm = .3micron, moly density = 10.2 g/cm^3 = 10.2e-9 mg/micron^3
    molyLoss = areaLoss*.3*10.2*(10**-6) #moly loss in micrograms
    molyLoss = round(molyLoss,2)
   	
    # Convert area values to mm^2 instead of microns^2, and round:
    PtAreaBef = round(PtAreaBef*10**-6, 4)
    PtAreaAft = round(PtAreaAft*10**-6, 4)
    areaLoss = round(areaLoss*10**-6, 4)
    # put all results into return tuples
    ret = (PtAreaBef, PtAreaAft, areaLoss, molyLoss,pcntMolyLoss)
    #picts = (Ptbef, Ptaft)
    print "Point"+I
    return ret

def CALisolatePt (image, P=90,D=120,M=180,Pt=180):
    """
    This function filters and thresholds the image using regionalThres in order 
    to estimate the area of exposed Pt. Argument "image" must be ndarray.
    """
    poster = makePoster(image)
    # Threshold the image. This is a global threshold. There is probably a better one out there.
    isoPt = regionalThresh(image, poster,p=P,d=D,m=M,pt=Pt,gaussBlur=3)
    # Make the image into a boolean image
    isoPt = isoPt.astype(np.bool_)
    return isoPt


def CaldirtComp (before, after, P=8,D=28,M=55,Pt=60):
    """
    This is the dirt compare foils function. It takes the pathnamess of two images (before
    and after) and returns the number of dirt particles and area of dirt on each foil, and 
    the thresholded images used to calculate the amount of dirt.   
    """
    # Load images
    bf = cv2.imread(before, cv2.CV_LOAD_IMAGE_GRAYSCALE)
    af = cv2.imread(after, cv2.CV_LOAD_IMAGE_GRAYSCALE)
	
    # Dirt analysis
    bposter,aposter = makePoster(bf),makePoster(af)
    bthreshed,bmasked = regionalThresh(bf, bposter, p=P,d=D,m=M,pt=Pt,MaskEdges=True,GetMask=True)
    athreshed,amasked = regionalThresh(af, aposter, p=P,d=D,m=M,pt=Pt,MaskEdges=True,GetMask=True)
    numbf,areabf,labbf,sizesbf = amountDirt(bthreshed)
    numaf,areaaf,labaf,sizesaf = amountDirt(athreshed)

    bthreshed = (bmasked/255)*bthreshed
    athreshed = (amasked/255)*athreshed

    # put all results into return tuples
    ret = (numbf,numaf,areabf,areaaf)
    #picts = (bthreshed, athreshed)
    return ret
    
def convTime(t):
    mins = round(t/60,0)
    seconds = round(t%60,3)
    return str(mins)+"mins "+str(seconds)+"secs"