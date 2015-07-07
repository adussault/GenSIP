import numpy as np
import GenSIP.functions as fun
import GenSIP.bigscans.bigfoils as bf

from matplotlib import pyplot as plt
import os
from socket import gethostname
from time import time
import pickle

###################################################################################

###################################################################################

def threshPrec(panFolder,maskFolder, MoDirt="Mo",**kwargs):
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
    # Get kwargs:
    foilname=kwargs.get("foilname","40360_2")
    Quarter=kwargs.get("quarter","Q1")
    Domain=kwargs.get("domain",40)
    step= kwargs.get("step",2)
    RUN=kwargs.get("RUN",("P","D","M","hE","Pt"))
    Pstart=kwargs.get("Pstart",130)
    Dstart=kwargs.get("Dstart",160)
    Mstart=kwargs.get("Mstart",190)
    hEstart=kwargs.get("hEstart",230)
    Ptstart=kwargs.get("Ptstart",240)
    res=kwargs.get("res",4)

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
        
    if "hE" in RUN:
        hE = np.linspace(hEstart,hEstart+Domain, numpoints).astype(np.uint16)
        hE[hE>255] = 255
    else:
        hE = np.zeros((numpoints)).astype(np.uint16)
        hE[:] = hEstart
        
    if "Pt" in RUN:
        Pt = np.linspace(Ptstart,Ptstart+Domain, numpoints).astype(np.uint16)
        Pt[Pt>255] = 255
    else:
        Pt = np.zeros((numpoints)).astype(np.uint16)
        Pt[:] = Ptstart
        

        
        
    # Generate plot subtitle:
    subtitle = "varies "
    Startstring = "from: "
    starts = {"P":Pstart,"D":Dstart,"M":Mstart, "hE":hEstart,"Pt":Ptstart}
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
    filetypes = ['.tif', '.jpg', '.jpeg','.tiff','.png']
    # Make list of files in the corresponding Before and After folders
    panSubs = os.listdir(panFolder)
    maskSubs = os.listdir(maskFolder)
    
    # Make output folder if does not exist
    outFolder = 'Output/Output_'+foilname+"_"+Quarter
    if not os.path.exists(outFolder):
        os.makedirs(outFolder)
    if not os.path.exists(outFolder+'/Plots/'+'Threshold_'+MoDirt):
        os.makedirs(outFolder+'/Plots/'+'Threshold_'+MoDirt)
    
    
    allowableMo = ("Mo","Moly","moly","molybdenum","Molybdenum","M","m")
    allowableDirt = ("Dirt","dirt","D","d")
    
    # Create a Dictionary to store the data:
    datestring = fun.getDateString()
    datestring.replace(':','.')
    host = os.path.splitext(gethostname())[0]
    Version = fun.getGenSIPVersion()
    
    Data = {
    "Date":datestring,"Host":host, 
    "Version":Version,"Analysis":MoDirt,"Domain":x,
    "SETTINGS":{"RUN":RUN,"Pstart":Pstart,"Dstart":Dstart,"Mstart":Mstart,\
    "hEstart":hEstart, "Ptstart":Ptstart},
    "DATA":{}
    }
    Data["Domain"] = x
    
    
 
                
    y1 = np.zeros((numpoints))
    y2 = np.zeros((numpoints))
    y3 = np.zeros((numpoints))
    
    if MoDirt in allowableMo:
        for i in range(numpoints):
            t1=time()
            ret = BIGanalyzeMoly(panSubs,panFolder, maskFolder, P=P[i],D=D[i],\
            M=M[i],HE=hE[i],Pt=Pt[i],I=str(i),MoDirt=MoDirt,res=res)
            totFoilArea = ret[0]
            totPtArea = ret[1]
            #PerPtLoss = (Ptbef-Ptaft)/Ptbef
            PercPt = ret[2]

            y1[i] = PercPt
            y2[i] = totFoilArea
            y3[i] = totPtArea
            t2 = time() - t1
            convTime(t2)
            
        
        # MAKE PLOT 
        
        plt.figure()
        plt.suptitle(foilname+" "+MoDirt+" Analysis")
        plt.title(subtitle,fontsize=14)
        plt.plot(x, y1, "-ro")
        plt.xlabel("Threshold Values")
        plt.ylabel("% Exposed Pt")
        #Set Max and min of the y axis:
        minY = min(np.min(y1),np.min(y2)) - 2
        maxY = max(np.max(y1),np.max(y2)) + 5
        if minY >= 0:
            minY = 0
        plt.axis([x[0]-2,x[numpoints-1]+2,minY, maxY])
        plt.savefig(outFolder+'/Plots/'+'Threshold_'+MoDirt+"_"+datestring+"/"+foilname+" Mo Threshold plot.png")
        print "Percent Pt plot for "+foilname+" finished."
        plt.close()
        
        plt.figure()
        plt.suptitle(foilname+" "+MoDirt+" Analysis")
        plt.title(subtitle,fontsize=14)
        plt.plot(x, y2, "-ro")
        plt.plot(x, y3, "-bo")
        plt.legend(["Foil Area","Pt Area"],loc=0)
        plt.xlabel("Threshold Values")
        plt.ylabel("Area (cm^2)")
        #Set Max and min of the y axis:
        minY = min(np.min(y2),np.min(y3)) - 5
        maxY = max(np.max(y2),np.max(y3)) + 5
        if minY >= 0:
            minY = 0
        plt.axis([x[0]-2,x[numpoints-1]+2,minY, maxY])
        plt.savefig(outFolder+'/Plots/'+'Threshold_'+MoDirt+"_"+datestring+"/"+foilname+" Pt area plot.png")
        print "Pt area plot for "+foilname+" finished."
        plt.close()
        
        plt.figure()
        plt.suptitle(foilname+" "+MoDirt+" Derivative")
        plt.title(subtitle,fontsize=14)
        grad1 = np.gradient(y1, step)
        grad2 = np.gradient(y2, step)
        plt.plot(x, grad1, "-ro", label = "Deriv Dirt loss by count")
        plt.plot(x, grad2, "-bo", label="Deriv Dirt loss by area")
        plt.legend(["Deriv mass Mo Lost","Deriv % Mo Loss"],loc=0)
        plt.xlabel("Threshold Values")
        plt.ylabel("dy/dx")
        #Set Max and min of the y axis:
        minY = min(np.min(grad1),np.min(grad2)) - 5
        maxY = max(np.max(grad1),np.max(grad2)) + 5
        if minY >= 0:
            minY = 0
        plt.axis([x[0]-2,x[numpoints-1]+2,minY, maxY])
        plt.savefig(outFolder+'/Plots/'+'Threshold_'+MoDirt+"_"+datestring+"/"+foilname+" Mo Gradient plot.png")
        plt.close()
        
        Data["DATA"][foilname] = {"PercPt":y1,"AreaFoil":y2,"AreaPt":y3,"percGrad":grad1,"AreaGrad":grad2}
        
        print "Moly gradient plot for "+foilname+" finished."
                    
    elif MoDirt in allowableDirt: ###<----------NOT YET USABLE
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
        print "Dirt plot for "+foilname+" finished."
        
        plt.figure()
        plt.suptitle(foilname+" "+MoDirt+" Analysis")
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
        plt.savefig(outFolder++'/Plots/'+'Threshold_'+MoDirt+"_"+datestring+"/"+foilname+" Dirt Threshold plot.png")
        plt.close()
        
        plt.figure()
        plt.suptitle(foilname+" "+MoDirt+" Derivative")
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
        plt.savefig(outFolder+'/Plots/'+'Threshold_'+MoDirt+"_"+datestring+"/"+foilname+" Dirt Gradient plot.png")
        plt.close()
        
        plt.figure()
        plt.suptitle(foilname+" "+MoDirt+" Analysis")
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
        plt.savefig(outFolder+'/Plots/'+'Threshold_'+MoDirt+"_"+datestring+"/"+foilname+" Dirt Threshold plot.png")
        plt.close()
        
    else:
        raise Exception("\n MoDirt must be one of the following: \n"+\
        str(allowableMo)+"\n - - - - - or - - - - - \n"+str(allowableDirt))

    # Get the run time of the function and add it to the dictionary "Data":
    t2 = time()-t1
    Data["TimeEllapse"] = convTime(t2)
    # Serialize the "Data" dictionary and store it in a pickle file with the other images. 
    with open(outFolder+'/Plots/'+'Threshold_Data_'+MoDirt+str(Domain)+"_"+datestring,'wb') as f:
        pickle.dump(Data,f)

###################################################################################

###################################################################################

def BIGanalyzeMoly(panSubs,panFolder,maskFolder, P=150,D=180,M=210,HE=253,Pt=240,I='',MoDirt='Mo',res=4):
    panSubs=os.listdir(panFolder)
    panSubs[:] = [img for img in panSubs if img.startswith("sub_")]
    panSubs[:] = [img for img in panSubs if img.endswith(".tif")]
    totPtArea = 0
    totFoilArea = 0
    for sub in panSubs:
        root,ext = os.path.splitext(sub)
        
        # Create the threshholded image
        subImage = fun.loadImg(panFolder+'/'+sub,0)
        subMask = fun.loadImg(maskFolder+'/'+sub,0)
        proc = bf.bigPostPreProc(subImage)
        poster = bf.bigPosterfy(proc)
        threshed = bf.bigRegionalThresh(subImage,poster,Mask=subMask,\
        p=P,d=D,m=M,hE=HE,pt=Pt,MoDirt=MoDirt)
        threshed = threshed.astype(np.bool_)
        
        # Get amount of exposed Platinum:
        
        PixPt = np.sum(threshed)
        AreaPt = round(PixPt*res*10**-6, 4)
        PixFoil = np.sum(subMask.astype(np.bool_))
        AreaFoil = round(PixFoil*res*10**-6, 4)
        
        if AreaFoil == 0:
            PercPt = 0
        else:
            PercPt = round(float(AreaPt)/float(AreaFoil)*100,2)
        
        # Make output image
        threshed=threshed.astype(np.uint8)*255
        
        totPtArea = totPtArea + AreaPt
        totFoilArea = totFoilArea + AreaFoil
        del(threshed)
        
    totPtArea = round(float(totPtArea)/100,2)
    totFoilArea = round(float(totFoilArea)/100,2)
    if totFoilArea == 0:
        PercPt = np.nan
    else:
        PercPt = round(100*totPtArea/totFoilArea,2)
    print "Point"+I
    return (totFoilArea,totPtArea,PercPt)
    
###################################################################################

###################################################################################

def convTime(t):
    mins = round(t/60,0)
    seconds = round(t%60,3)
    print str(mins)+"mins "+str(seconds)+"secs"
    
#Q1 = fun.loadImg("InputPicts/FoilScans/Q1/panorama.tif",0)
