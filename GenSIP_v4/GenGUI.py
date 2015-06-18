import Tkinter as Tk
import cv2
import matplotlib.pyplot as plt
import mahotas as mh
from scipy import misc
import GenSIP_v4.functions as fun
import GenSIP_v4.Kuwahara as Kuwahara
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
import numpy as np
#from GenSIP_v4.loadTestImages import *


def GUIfy(image):
    """
    This method takes an image of a foil (a numpy ndarray, not a pathname),and starts 
    up a GUI that allows the user to adjust the settings on various functions in 
    order to see their effects on the foil. The most developed of the sub-GUIs is 
    the regionalThresh GUI. 
    """
    def makeOdd(n):
        #global past
        n = int(n)
        if not n%2:
            Tk.blocksize.set(n+1)
    
    def win_destroy(root):
        root.destroy()
    def r1Destroy():
        root1.destroy()
    ########################################################
    ################       SETUP       #####################
    ########################################################
    def setupGUI():
        win_destroy(root1)
        global root2
        global f
        global a
        global canvas
        global toolbar
        root2=Tk()
        f = Figure()
        a = f.add_subplot(111)
        a.imshow(image,"gray")
        a.set_title("40393,0230 before cleaning")
        
        canvas = FigureCanvasTkAgg(f, master=root2)
        canvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH,expand=1)
        canvas._tkcanvas.grid(row=0,column=0,columnspan=3,sticky=Tk.SW)
        #Navigation Toolbar automatically calls pack, so it requires a frame,
        # which can then be put in a grid:
        toolbar_frame = Tk.Frame(root2)
        toolbar_frame.grid(row=1,column=0,columnspan=3)
        toolbar = NavigationToolbar2TkAgg(canvas, toolbar_frame)   
        toolbar.update()
        
    def resetGUI():
        root2.destroy()
        GUIfy(image)
        
    def showOG():
            a.imshow(image,"gray")
            a.figure.canvas.draw()
            
    ####________________________________________________________####
    ####
    ####________________________________________________________####
    
    
    #############################################################
    ########                   GUICONTS                    ######
    #############################################################
    def GUIconts():
        
        def FindContours():
            t1 = Thresh1.get()
            t2 = Thresh2.get()
            canny = cv2.Canny(image, t1,t2)
            a.imshow(canny,"gray")
            a.figure.canvas.draw()
            
        setupGUI()
        
        Thresh1=Tk.Scale(root2,from_=0, to=255,resolution=1,orient=Tk.HORIZONTAL,length=255,\
        label='Threshold 1')
        Thresh1.set(55)
        Thresh1.grid(row=2,column=0,columnspan=2)
        
        Thresh2 = Tk.Scale(root2, from_=0,to=255, resolution=1,orient=Tk.HORIZONTAL,length=255,\
        label='Threshold 2')
        Thresh2.set(200)
        Thresh2.grid(row=3,column=0)
        
        Update = Tk.Button(root2, text="Update", command=FindContours)
        Update.grid(row=4,column=0)
        
        CLOSE = Tk.Button(root2, text="Main Menu", command=resetGUI)
        CLOSE.grid(row=4,column=1)
        root2.mainloop()
        
    #############################################################
    ########                    GUIADPTRSH                 ######
    #############################################################
    def GUIadptrsh():
        def AdaptThresh():
            maVal = 255#maxVal.get()
            bsize = blocksize.get()
            c = C.get()
            adtv = cv2.adaptiveThreshold(image,maVal,cv2.ADAPTIVE_THRESH_GAUSSIAN_C,cv2.THRESH_BINARY,\
            bsize,c)
            a.imshow(adtv,"gray")
            a.figure.canvas.draw()
            
        setupGUI()
        blocksize=Tk.Scale(root2,from_=2, to=1001,resolution=1,orient=Tk.HORIZONTAL,length=500,\
        label='Block Size',command=makeOdd)
        blocksize.set(55)
        #blocksize.pack(side=LEFT, after=maxVal)
        blocksize.grid(row=2,column=0,columnspan=2)
        
        C = Tk.Scale(root2, from_=0,to=255, label="C",orient=Tk.HORIZONTAL,\
        length=255)
        C.set(2)
        #C.pack(side=RIGHT)
        C.grid(row=3,column=0)
        
        Update = Tk.Button(root2, text="Update", command=AdaptThresh)
        Update.grid(row=4,column=0)
        
        CLOSE = Tk.Button(root2, text="Main Menu", command=resetGUI)
        CLOSE.grid(row=4,column=1)
        
        root2.mainloop()
    #############################################################
    ########                   GUIBILAT                    ######
    #############################################################
    def GUIbilat():
        def Bilat():
            d = D.get()
            sigColor=sigC.get()
            sigSpace=sigS.get()
            bilat = cv2.bilateralFilter(image, d, sigColor, sigSpace)
            a.imshow(bilat,"gray")
            a.figure.canvas.draw()
            
        setupGUI()
        D=Tk.Scale(root2,from_=-10, to=500,resolution=1,orient=Tk.HORIZONTAL,length=500,\
        label='Distance')
        D.set(55)
        D.grid(row=2,column=0,columnspan=2)
        
        sigC = Tk.Scale(root2, from_=0,to=255, label="SigColor",orient=Tk.HORIZONTAL,\
        length=255)
        sigC.set(10)
        sigC.grid(row=3,column=0)
        
        sigS = Tk.Scale(root2, from_=0, to=500, label="SigSpace",orient=Tk.HORIZONTAL,\
        length=250)
        sigS.set(20)
        sigS.grid(row=3,column=1)
        
        Update = Tk.Button(root2, text="Update", command=Bilat)
        Update.grid(row=4,column=0)
        CLOSE = Tk.Button(root2, text="Main Menu", command=resetGUI)
        CLOSE.grid(row=4,column=1)
        root2.mainloop()
    #############################################################
    ########                   GUIPOST                     ######
    #############################################################
    def GUIpost():
        def pfy():
            k = K.get()
            p = posterfy(image,k)
            a.imshow(p,"gray")
            a.figure.canvas.draw()
        def cfy():
            t1 = Thresh1.get()
            t2 = Thresh2.get()
            k = K.get()
            p = posterfy(image,k)
            canny = cv2.Canny(p, t1,t2)
            a.imshow(canny,"gray")
            a.figure.canvas.draw()
        def ffy():
            t1 = Thresh1.get()
            t2 = Thresh2.get()
            k = K.get()
            p = posterfy(image,k)
            ncs=NumConts.get()
            canny = cv2.Canny(p, t1,t2)
            conts,hier = cv2.findContours(canny.copy(),cv2.RETR_TREE,cv2.CHAIN_APPROX_NONE)
            h = fillAllConts(canny,conts,ncs)
            a.imshow(h)
            a.figure.canvas.draw()
        def posterfy(ar,k_size=5):
            nar=ar.copy()
            blk=(nar<=4)
            pleat=(5<=nar)&(nar<=75)
            darkMo=(76<=nar)&(nar<=99)
            Mo=(100<=nar)&(nar<=199)
            Pt=(200<=nar)
            nar[blk]=0
            nar[pleat]=50
            nar[darkMo]=85
            nar[Mo]=150
            nar[Pt]=255
            nar = fun.maskEdge(nar,thickness=30)
            #Do a morphological opening step to eliminate the rough edges
            # I prefer opening over closing because it is more important to catch
            # all of the dark areas, as the Pt areas will have pretty sraightforward
            # threshold results
            kernel = np.ones((k_size,k_size))
            nar = cv2.morphologyEx(nar, cv2.MORPH_OPEN, kernel)
            return nar
        def fillAllConts(img,conts,numconts):
            img = cv2.cvtColor(img,cv2.COLOR_GRAY2RGB)
            colors = []
            if F:
                #fill all contours if fill box is checked
                opt=-1
            else:
                #draw contour perimeter with thickness of 5 otherwise
                opt=5
            if numconts > len(conts):
                nc = range(len(conts))
            else:
                nc=range(numconts)
            for x in nc:
                colors.append((x*7,255-x*7,0))
            for i in nc:
                cv2.drawContours(img,conts,i,colors[i],opt)
            return img
        
        setupGUI()
        K=Tk.Scale(root2,from_=3, to=20,resolution=1,orient=Tk.HORIZONTAL,length=80,\
        label='kernel size')
        K.set(6)
        K.grid(row=2,column=0)
        
        Thresh1=Tk.Scale(root2,from_=0, to=255,resolution=1,orient=Tk.HORIZONTAL,length=255,\
        label='Threshold 1')
        Thresh1.set(55)
        Thresh1.grid(row=2,column=1,columnspan=1)
        
        Thresh2 = Tk.Scale(root2, from_=0,to=255, resolution=1,orient=Tk.HORIZONTAL,length=255,\
        label='Threshold 2')
        Thresh2.set(200)
        Thresh2.grid(row=3,column=0)
        
        NumConts = Tk.Scale(root2, from_=0,to=50, resolution=1,orient=Tk.HORIZONTAL,length=100,\
        label='Number of Contours to examine')
        NumConts.set(33)
        NumConts.grid(row=4,column=0)
        F=Tk.BooleanVar()
        Fill = Tk.CheckTk.Button(root2,text="Fill Contours",variable=F)

        Fill.grid(row=3,column=1)
        
        pofy = Tk.Button(root2, text="Posterfy", command=pfy)
        pofy.grid(row=5,column=1)
        
        contfy = Tk.Button(root2, text="Get Contours", command=cfy)
        contfy.grid(row=5,column=0)
        
        fillfy = Tk.Button(root2, text="Paint Contours", command=ffy)
        fillfy.grid(row=4,column=1,columnspan=1)
        
        ORIG = Tk.Button(root2, text="Show Original", command=showOG)
        ORIG.grid(row=6,column=0)
        
        CLOSE = Tk.Button(root2, text="Main Menu", command=resetGUI)
        CLOSE.grid(row=6,column=1)
        root2.mainloop()
    ##############################################################
    ########                   GUITHRESH                    ######
    ##############################################################
    def GUIthresh():
        def Thresh():
            maVal = MAX.get()
            thresh = THRESH.get()
            V = long(v.get())
            gauss = BLUR.get()
            imageF = cv2.GaussianBlur(image, (gauss,gauss), 0)
            ret,thr = cv2.threshold(imageF,thresh,maVal,V)
            a.imshow(thr,"gray")
            a.figure.canvas.draw()
        def cann():
            maVal = MAX.get()
            thresh = THRESH.get()
            V = long(v.get())
            ret,thr = cv2.threshold(image,thresh,maVal,V)
            can = cv2.Canny(thr,thresh,maVal)
            a.imshow(can,"gray")
            a.figure.canvas.draw()
        
        def makeGaussOdd(n):
        #global past
            n = int(n)
            if not n%2:
                BLUR.set(n+1)
                
        setupGUI()
        
        THRESH = Tk.Scale(root2, from_=0,to=255, label="Threshold",orient=Tk.HORIZONTAL,\
        length=255)
        THRESH.set(2)
        THRESH.grid(row=2,column=0,columnspan=2)
        
        MAX = Tk.Scale(root2, from_=0,to=255, label="Maximum Value",orient=Tk.HORIZONTAL,\
        length=255)
        MAX.set(255)
        MAX.grid(row=3,column=0,columnspan=2)
        
        BLUR = Tk.Scale(root2, from_=0,to=15, label="Gaussian Blur",orient=Tk.HORIZONTAL,\
        length=60,command=makeGaussOdd)
        BLUR.set(3)
        BLUR.grid(row=2,column=2,columnspan=1)
        
        v = Tk.IntVar()
        Tk.RadioTk.Button(root2, text="Binary", variable=v, value=0).grid(row=4,column=0)
        Tk.RadioTk.Button(root2, text="Binary Inverse", variable=v, value=1).grid(row=4,column=1)
        Tk.RadioTk.Button(root2, text="Truncated", variable=v, value=2).grid(row=5,column=0)
        Tk.RadioTk.Button(root2, text="To Zero", variable=v, value=3).grid(row=5,column=1)
        Tk.RadioTk.Button(root2, text="To Zero Inverse", variable=v, value=4).grid(row=6,column=0)
        Tk.RadioTk.Button(root2, text="Mask", variable=v, value=7).grid(row=6,column=1)
        Tk.RadioTk.Button(root2, text="Otsu", variable=v, value=8).grid(row=7,column=0)
        
        
        Update = Tk.Button(root2, text="Threshold", command=Thresh)
        Update.grid(row=3,column=2)
        
        CAN = Tk.Button(root2, text="Find Edges", command=cann)
        CAN.grid(row=4,column=2)
        
        ORIG = Tk.Button(root2, text="Show Original", command=showOG)
        ORIG.grid(row=5,column=2)
        
        CLOSE = Tk.Button(root2, text="Main Menu", command=resetGUI)
        CLOSE.grid(row=7,column=2)
        
        root2.mainloop()
        
    #############################################################
    ########               REGIONAL THRESH                 ######
    #############################################################
    def GUIreg():
    # This GUI has several functions plugged into it: posterfy, regionalThresh,
    #  Threshold, and pre-Poster.
        def pfy():
            k = K.get()
            kuw = KUW.get()
            g1=G1.get()
            g2=G2.get()
            global poster
            poster = fun.makePoster (image,k,kuw,g1,g2,.1)
            a.imshow(poster,"gray")
            a.figure.canvas.draw()
            
        def prePost():
            KuSize = KUW.get()
            Gaus1=G1.get()
            Gaus2=G2.get()
            rsz = misc.imresize(image,.1,interp='bicubic')
            gr = cv2.GaussianBlur(rsz, (Gaus1,Gaus1),0)
            kgr = Kuwahara(gr,KuSize)
            rkgr = misc.imresize(kgr,(image.shape),interp='bicubic')
            grkgr = cv2.GaussianBlur(rkgr, (Gaus2,Gaus2),0)
            a.imshow(grkgr,"gray")
            a.figure.canvas.draw()
            
        def UpdatePoster():
            k = K.get()
            kuw = KUW.get()
            g1=G1.get()
            g2=G2.get()
            global poster
            poster = fun.makePoster(image,k,kuw,g1,g2,.1)
            
        def calcDirt():
            p=PLEAT.get()
            d=DARKMO.get()
            m=MO.get()
            pt=PT.get()
            gauss=BLUR.get()
            thrshType=long(v.get())
            thisposter = fun.makePoster(image)
            thisthreshed = fun.regionalThresh(image,thisposter,p,d,m,pt,gauss,thrshType,MaskEdges=True)
            inv = cv2.bitwise_not(thisthreshed)
            inv = inv.astype(np.uint8)
            labelledFoil,numDirt = mh.label(inv)
            areaDirt = np.sum(thisthreshed/255)
            print "Number of Dirt Particles: %s" % numDirt
            print "Area of Dirt Particles: %s" % areaDirt
            print "Inv %s " % inv.dtype
            print inv[0,0]
            print "labelled: %s" %labelledFoil.dtype
            print "Threshed: %s"%thisthreshed.dtype
            print thisthreshed[0,0] 
            a.imshow(labelledFoil,"jet")
            a.figure.canvas.draw()
            
        def calcPt():
            UpdatePoster()
            p=PLEAT.get()
            d=DARKMO.get()
            m=MO.get()
            pt=PT.get()
            gauss=BLUR.get()
            thrshType=long(v.get())
            global threshed
            threshed = fun.regionalThresh(image,poster,p,d,m,pt,gauss,thrshType)
            a.imshow(threshed,"gray")
            a.figure.canvas.draw()
            binthresh = threshed.astype(np.bool_)
            PtArea = fun.calcPtArea(binthresh)
            print "Area of Pt (pixels): %s" % PtArea
            
        def fT():
            UpdatePoster()
            p=PLEAT.get()
            d=DARKMO.get()
            m=MO.get()
            pt=PT.get()
            gauss=BLUR.get()
            thrshType=long(v.get())
            global threshed
            threshed = regionalThresh(image,poster,p,d,m,pt,gauss,thrshType,True)
            a.imshow(threshed,"gray")
            a.figure.canvas.draw()
        
        def Thresh():
            maVal = 255
            thresh = MO.get()
            V = long(v.get())
            gauss = BLUR.get()
            imageF = cv2.GaussianBlur(image, (gauss,gauss), 0)
            ret,thr = cv2.threshold(imageF,thresh,maVal,V)
            a.imshow(thr,"gray")
            a.figure.canvas.draw()
            
        def cann():
            thr=MO.get()
            UpdatePoster()
            can = cv2.Canny(thr,50,255)
            a.imshow(can,"gray")
            a.figure.canvas.draw()
        
        def makeGaussOdd(n):
        #global past
            n = int(n)
            if not n%2:
                BLUR.set(n+1)
        def makeG1Odd(n):
        #global past
            n = int(n)
            if not n%2:
                G1.set(n+1)
        def makeG2Odd(n):
        #global past
            n = int(n)
            if not n%2:
                G2.set(n+1)
        def fixKuwahara(n):
        #global past
            n = int(n)
            if n%4!=1:
                KUW.set(n-n%4+1)
                
        setupGUI()

        ###################################
        ####        VARIABLES          ####
        ###################################
        
            ### Threshold variables ###
        PLEAT = Tk.Scale(root2, from_=0,to=255, label="Pleat",orient=Tk.HORIZONTAL,\
        length=255)
        PLEAT.set(5)
        PLEAT.grid(row=2,column=0,columnspan=2)
        
        DARKMO = Tk.Scale(root2, from_=0,to=255, label="Dark Moly",orient=Tk.HORIZONTAL,\
        length=255)
        DARKMO.set(25)
        DARKMO.grid(row=3,column=0,columnspan=2)
        
        MO = Tk.Scale(root2, from_=0,to=255, label="Moly",orient=Tk.HORIZONTAL,\
        length=255)
        MO.set(55)
        MO.grid(row=4,column=0,columnspan=2)
        
        PT = Tk.Scale(root2, from_=0,to=255, label="Platinum",orient=Tk.HORIZONTAL,\
        length=255)
        PT.set(60)
        PT.grid(row=5,column=0,columnspan=2)
        
        BLUR = Tk.Scale(root2, from_=0,to=15, label="Gaussian Blur",orient=Tk.HORIZONTAL,\
        length=60,command=makeGaussOdd)
        BLUR.set(3)
        BLUR.grid(row=2,column=2,columnspan=1)
            ### Posterfy variables ###
        K=Tk.Scale(root2,from_=3, to=20,resolution=1,orient=Tk.HORIZONTAL,length=80,\
        label='Kernel size')
        K.set(6)
        K.grid(row=3,column=2)
            ### makePoster Variables ###
        KUW=Tk.Scale(root2,from_=3, to=21,resolution=1,orient=Tk.HORIZONTAL,length=84,\
        label='Kuwahara size',command=fixKuwahara)
        KUW.set(9)
        KUW.grid(row=2,column=4)
        
        G1=Tk.Scale(root2, from_=0,to=25, label="1st Gaussian Blur",orient=Tk.HORIZONTAL,\
        length=60,command=makeG1Odd)
        G1.set(3)
        G1.grid(row=3,column=4)
        
        G2=Tk.Scale(root2, from_=0,to=25, label="2nd Gaussian Blur",orient=Tk.HORIZONTAL,\
        length=60,command=makeG2Odd)
        G2.set(11)
        G2.grid(row=4,column=4)
        
        
        ####### Tk.RadioTk.ButtonS #######
        
        v = Tk.IntVar()
        Tk.RadioTk.Button(root2, text="Binary", variable=v, value=0).grid(row=1,column=3)
        Tk.RadioTk.Button(root2, text="Binary Inverse", variable=v, value=1).grid(row=2,column=3)
        Tk.RadioTk.Button(root2, text="Truncated", variable=v, value=2).grid(row=3,column=3)
        Tk.RadioTk.Button(root2, text="To Zero", variable=v, value=3).grid(row=4,column=3)
        Tk.RadioTk.Button(root2, text="To Zero Inverse", variable=v, value=4).grid(row=5,column=3)
        Tk.RadioTk.Button(root2, text="Otsu", variable=v, value=8).grid(row=6,column=3)
        
        ###############################
        #######     Tk.ButtonS     #######
        ###############################
        
        pofy = Tk.Button(root2, text="Posterfy", command=pfy)
        pofy.grid(row=4,column=2)
        
        FT = Tk.Button(root2, text="Regional Thresh", command=fT)
        FT.grid(row=5,column=2)
        
        ORIG = Tk.Button(root2, text="Show Original", command=showOG)
        ORIG.grid(row=6,column=2)
        
        DIRTcal = Tk.Button(root2, text="Print Amount Dirt", command=calcDirt)
        DIRTcal.grid(row=6,column=4) 
        
        PTcal = Tk.Button(root2, text="Print Pt Area", command=calcPt)
        PTcal.grid(row=7,column=4) 
                  
        PROC = Tk.Button(root2, text="Show Pre-Poster", command=prePost)
        PROC.grid(row=5,column=4)      
        
        CAN = Tk.Button(root2, text="Find Edges", command=cann)
        CAN.grid(row=7,column=2)
        
        TH = Tk.Button(root2, text="Threshold", command=Thresh)
        TH.grid(row=7,column=3)
        
        CLOSE = Tk.Button(root2, text="Main Menu", command=resetGUI)
        CLOSE.grid(row=7,column=0)

        
        root2.mainloop()
        
    # Startup window
    root1 = Tk()
    FC = Tk.Button(root1, text="FindContours", command=GUIconts)
    FC.grid(row=0,column=0)
    AT = Tk.Button(root1, text="Adaptive Threshold", command=GUIadptrsh)
    AT.grid(row=1, column=0)
    BIL = Tk.Button(root1, text="Bilateral Filter", command=GUIbilat)
    BIL.grid(row=2, column=0)
    PO = Tk.Button(root1,text="Region Selection", command=GUIpost)
    PO.grid(row=3,column=0)
    TSH = Tk.Button(root1,text="Thresholding", command=GUIthresh)
    TSH.grid(row=4,column=0)
    REG = Tk.Button(root1,text="Regional Thresh", command=GUIreg)
    REG.grid(row=5,column=0)
    CLOZ = Tk.Button(root1, text="Cancel", command=r1Destroy)
    CLOZ.grid(row=6,column=0)
    
    root1.mainloop()



'''
kuwa230 = cv2.imread("InputPicts/Before/Before_TEST/Kuwahara230.tif",cv2.CV_LOAD_IMAGE_GRAYSCALE)    
kuwa231 = cv2.imread("InputPicts/Before/Before_TEST/Kuwahara231.tif",cv2.CV_LOAD_IMAGE_GRAYSCALE)
kuwa232 = cv2.imread("InputPicts/Before/Before_TEST/Kuwahara232.tif",cv2.CV_LOAD_IMAGE_GRAYSCALE)

rekuwa231 = cv2.imread("InputPicts/Before/Before_TEST/0231 Resize Kuwahara.tif",cv2.CV_LOAD_IMAGE_GRAYSCALE)
Prekuwa231 = cv2.imread("InputPicts/Before/Before_TEST/0231 Resize Kuwahara.png",cv2.CV_LOAD_IMAGE_GRAYSCALE)
'''
