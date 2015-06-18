#This is the module that defines the functions used in Genesis SEM Image Processing (GenSIP)
#===============================================================
#========================###IMPORTS###==========================
#===============================================================
#import relevant modules
import cv2
import numpy as np
import matplotlib.pyplot as plt
import mahotas as mh
import csv
import os
from scipy import misc
from GenSIP.Kuwahara import Kuwahara


#===============================================================
#=======================###FUNCTIONS###=========================
#===============================================================
#_______________________________________________________________
#####################--MASTER FUNCTIONS--#######################
def analyzeMoly (sss, res=1):
# This is the master function for comparing the exposed Pt on foil images.
#  It takes in a 'Sample Set String' ('sss') and runs the comparison functions 
#  on all of the foil pictures in the before and after folders of that sample set. 
#  res is the resolution of the image, as micron per pixel. default set to one.
#  The Sample Set String is a short identifier for whichever set of SEM scans you are running.
    # Declare a list of acceptable file types: tif and jpg
	filetypes = ['.tif', '.jpg', '.jpeg','.tiff']
	# Make list of files in the corresponding Before and After folders
	befpics = sorted(os.listdir('InputPicts/Before/Before_'+sss))
	aftpics = sorted(os.listdir('InputPicts/After/After_'+sss))
	
	# Make output folder if does not exist
	if not os.path.exists('Output/Output_'+sss):
		os.makedirs('Output/Output_'+sss)
	if not os.path.exists('Output/Output_'+sss+'/PtMaps'):
		os.makedirs('Output/Output_'+sss+'/PtMaps')
	
	# Make dirt and Pt output csv files and corresponding writers
	Ptcsv = open('Output/Output_'+sss+'/Pt_output_'+sss+'.csv','w+b')
	Ptwriter = csv.writer(Ptcsv)
	# Make headers:
        Ptwriter.writerow(['Foil #','Pt Area Before (mm^2)', 'Pt Area After (mm^2)', \
	'Area of Mo Loss (mm^2)', 'Approx Mo Loss (micrograms)','% Mo lost'])
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
			print "Not a picture."
			
			
			
def analyzeDirt(sss):
# This is the master script that takes in a Sample Set String and runs the comparison 
#  functions on all of the foil pictures in the before and after folders of that sample set. 
#  The Sample Set String is a short identifier for whichever set of SEM scans you are running.
	# Declare a list of acceptable file types: tif and jpg
	filetypes = ['.tif', '.jpg', '.jpeg','.tiff']
	# Make list of files in the corresponding Before and After folders
	befpics = sorted(os.listdir('InputPicts/Before/Before_'+sss))
	aftpics = sorted(os.listdir('InputPicts/After/After_'+sss))
	
	# Make output folder if does not exist
	if not os.path.exists('Output/Output_'+sss):
		os.makedirs('Output/Output_'+sss)
	if not os.path.exists('Output/Output_'+sss+'/DirtMaps'):
		os.makedirs('Output/Output_'+sss+'/DirtMaps')

	# Make dirt and Pt output csv files and corresponding writers
	dirtcsv = open('Output/Output_'+sss+'/dirt_output_'+sss+'.csv','w+b')
	dirtwriter = csv.writer(dirtcsv)
	
	# Make headers:
	dirtwriter.writerow(['Foil #','Dirt Count Before', 'Dirt Count After', \
	'Dirt Area Before', 'Dirt Area After', 'Approx Percent Dirt Loss by Area'])

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
				dirtVals, dirtPicts =dirtComp(befpath,aftpath)
				
				(numbf,numaf,areabf,areaaf) = dirtVals
				(threshedbf, threshedaf) = dirtPicts
				# Calculate difference in area
				perDirtLoss = 100*(areaaf-areabf)/areabf
				
				# Write the output to a new line in the csv file
				dirtwriter.writerow([foilnum, numbf, numaf, areabf, areaaf, perDirtLoss])
				# Create the output images
				# outfoil = compareMap(picts)
				# Save outfoil to Output folder
				# cv2.imwrite('Output/Output_'+sss+'/'+foilnum+' compare'+exten, outfoil)
				# Save the threshed images for analysis
				cv2.imwrite('Output/Output_'+sss+'/DirtMaps/'+foilnum+' before_threshed.png',\
				 threshedbf, [cv2.cv.CV_IMWRITE_PNG_COMPRESSION,6])
				cv2.imwrite('Output/Output_'+sss+'/DirtMaps/'+foilnum+' after_threshed.png',\
				 threshedaf, [cv2.cv.CV_IMWRITE_PNG_COMPRESSION,6])
			else: print "No after image for foil "+foilnum
		else:
			print "Not a picture."

def dirtComp (before, after):
# This is the dirt compare foils function. It takes the pathnamess of two images (before
#  and after) and returns the number of dirt particles and area of dirt on each foil, and 
#  the thresholded images used to calculate the amount of dirt.   
	# Load images
	bf = cv2.imread(before, cv2.CV_LOAD_IMAGE_GRAYSCALE)
	af = cv2.imread(after, cv2.CV_LOAD_IMAGE_GRAYSCALE)
	
	# Dirt analysis
	bposter,aposter = makePoster(bf),makePoster(af)
	bthreshed,bmasked = regionalThresh(bf, bposter, MaskEdges=True,GetMask=True)
	athreshed,amasked = regionalThresh(af, aposter, MaskEdges=True,GetMask=True)
	numbf,areabf,labbf,sizesbf = amountDirt(bthreshed)
	numaf,areaaf,labaf,sizesaf = amountDirt(athreshed)
        
        bthreshed = (bmasked/255)*bthreshed
        athreshed = (amasked/255)*athreshed
        
	# put all results into return tuples
	ret = (numbf,numaf,areabf,areaaf)
	picts = (bthreshed, athreshed)
	return ret, picts

def MoComp (before, after, res=1):
# This is the exposed platinum compare foils function. It takes the pathnamess of two images (before
#  and after) and the image resolution, in square microns per pixel. Default set to 1.
#  NOTE: This function assumes before and after images have the same resolution!
# MoComp returns two tuples: 
#       -The first tuple contains the area of exposed platinum before and after 
#        cleaning, the area of moly loss, the micrograms of moly lost, and the 
#        approximate % of the original molybdenum lost.
#       -The second tuple contains the binary before and after images of the
#        exposed platinum used in these calculations. 

	# Load images
	bf = cv2.imread(before, cv2.CV_LOAD_IMAGE_GRAYSCALE)
	af = cv2.imread(after, cv2.CV_LOAD_IMAGE_GRAYSCALE)

	# Generate binary thresholds for Platinum:
	Ptbef,Ptaft = isolatePt(bf),isolatePt(af)
	
	# Approximate the percent of molybdenum lost:
	# Will take the average 
	foilareabef = getFoilArea(bf,res)
	foilareaaft = getFoilArea(af,res)
	# Calculate area of exposed platinum, adjust for resolution:
	PtAreaBef = calcPtArea(Ptbef)*res
	PtAreaAft = calcPtArea(Ptaft)*res
	
	# Calculate difference in area of exposed Pt:
	# Normalize values to the area of the before image
	befaftRatio = foilareabef/foilareaaft
	areaLoss = PtAreaAft*befaftRatio-PtAreaBef
        MolyBef = foilareabef-PtAreaBef
        MolyAft = foilareaaft-PtAreaAft
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
	picts = (Ptbef, Ptaft)
	return ret, picts


#_____________________________________________________________
######################--DIRT COUNTING--#########################

def maskEdge(img, thickness = 80):
# This function masks off the outer edge of the foil, since the edge complicates dirt particle counting
# Parameters are:
#	img = nparray of image,
#	mask_thickness = the number of pixels you want to take of the edge of the foil
# 		Based on the image resolution, 1 micron ~= i pixel
	#Load image
	image = img.copy()
	# Apply a morphological closing to the image to remove small dark areas:
	kernel = np.ones((3,3))
        cimage = cv2.morphologyEx(image, cv2.MORPH_CLOSE, kernel)
	#Threshold and blur image for contour detection
	blur2 = cv2.blur(cimage, (100,100))
	ret,thresh2 = cv2.threshold(blur2,50,255,cv2.THRESH_BINARY)
	# Put a black border around the thresholded image to prevent edge problems
        thresh2 = cv2.copyMakeBorder(thresh2,5,5,5,5,cv2.BORDER_CONSTANT,value=0)
	#Find the contours of the foil
	can1 = cv2.Canny(thresh2, 50, 200)
	contours,hierarchy = cv2.findContours(can1.copy(),cv2.RETR_TREE,cv2.CHAIN_APPROX_NONE)

	#Isolate and draw the largest contour, which should be the outline of the foil
	#Sort contours by size
	contours = sorted(contours, key = cv2.contourArea, reverse = True)
	#Largest contour should be the outer edge of foil
	outercnt = contours[0]
	#Create corresponding lists of contour areas and perimeters. 
	#cntAreas = list()
	#cntPerim = list()
	#Black out the edges of the foil
	#Double the mask_thickness to get the thickness of the masking line. 
	mask_thickness = thickness*2
	cv2.drawContours(image, [outercnt], -1, (0,255,0), mask_thickness)
	cv2.drawContours(thresh2, [outercnt], -1, (0,255,0), mask_thickness)
	# remove border from thresh2
	thresh2 = thresh2[5:thresh2.shape[0]-5,5:thresh2.shape[1]-5]
	# Return botht the image and the foil area for the sake of future moly loss
	# approximations and for dirt counting.
	return image,thresh2

	
def amountDirt(img):
# Calculates the number of dirt particles and the area of the foil covered by dirt
# Takes a black and white image (black dirt on white foil)
# Returns the number of dirt particles, the area of dirt particles, and a labelled image
	#inverse the colors since mahotas.label only works on white on black
	inv = cv2.bitwise_not(img)
	#invDirt = cv2.bitwise_not(isoDirt(img,profile))
	labeledFoil,numDirt = mh.label(inv)
	# Don't count the foil in numDirt:
	numDirt -= 1
	#Calculate the area of the dirt using Findcontours
	sizes = mh.labeled.labeled_size(labeledFoil)
	# Sort sizes of particles by size in descending order:
	sizes = np.sort(sizes)[::-1]
	# Eliminate the foil from the "sizes" array:
	sizes = sizes[2:]
	# Total area of dirt is equal to the sum of the sizes.
	areaDirt = sum(sizes)
	return (numDirt, areaDirt, labeledFoil, sizes)
	
def makePoster(image,kern=6, KuSize=9,Gaus1=3,Gaus2=11,rsize=.1):
# This method takes the image of the foil and creates a smoothed Kuwahara image
#  used to make the poster for regional thresholding.
    rsz = misc.imresize(image,rsize,interp='bicubic')
    gr = cv2.GaussianBlur(rsz, (Gaus1,Gaus1),0)
    kgr = Kuwahara(gr,KuSize)
    rkgr = misc.imresize(kgr,(image.shape),interp='bicubic')
    grkgr = cv2.GaussianBlur(rkgr, (Gaus2,Gaus2),0)
    prkgr = posterfy(grkgr,kern)
    return prkgr

def posterfy(image,k_size=6):
# Takes a gray image and sets all values within a given range to a single value
# this way we can divide up regions into primarily Pt, primarily Mo, dirt, or black
    image_copy=image.copy()
    image_copy=image_copy.astype(np.uint8)
    blk=(image_copy<=4)
    pleat=(5<=image_copy)&(image_copy<=75)
    darkMo=(76<=image_copy)&(image_copy<=110)
    Mo=(111<=image_copy)&(image_copy<=199)
    Pt=(200<=image_copy)
    image_copy[blk]=0
    image_copy[pleat]=50
    image_copy[darkMo]=85
    image_copy[Mo]=150
    image_copy[Pt]=255
    #image_copy = maskEdge(image_copy,thickness=60)
    #Do a morphological opening step to eliminate the rough edges
    # I prefer opening over closing because it is more important to catch
    # all of the dark areas, as the Pt areas will have pretty sraightforward
    # threshold results
    kernel = np.ones((k_size,k_size))
    image_copy = cv2.morphologyEx(image_copy, cv2.MORPH_OPEN, kernel)
    return image_copy
    
def regionalThresh(ogimage,poster,p=8,d=28,m=55,pt=60,gaussBlur=3,threshType=0L,MaskEdges=0,GetMask=0):
# This is the main thresholding method for analyzing the amount of Pt and dirt on
#   the foils. It takes the posterized image and splits the image into regions based on
#   the gray level in these different regions, then it assigns different threshold parameters 
#   to the different regions, specified by the arguments p,d,m, and pt, which correspond to the
#   regions 'pleat','dark Molybdenum', 'Molybdenum,' and 'platinum.'
# regionalThresh returns the final thresholded image with the black area around the foil
#   cut out so that only the dirt appears. This allows mh.label to count the dirt and not get
#   thrown off by the foil outline. If the option "GetMask" is set to True, then regionalThresh
#   also returns the image of the outline of the foil and all regions that are black (<5).
    if poster.shape != ogimage.shape:
        raise Exception("The two arrays are not the same shape.")
        return
    Image = ogimage.astype(np.uint8)
    gPoster = poster.astype(np.uint8)
    threshedImage = np.zeros((Image.shape),dtype=np.uint8)

    # Create the images
    blk = gPoster.copy()
    pleat = gPoster.copy()
    darkMo = gPoster.copy()
    Mo = gPoster.copy()
    Pt = gPoster.copy()
    # isolate various sections
    blk[blk!=0]=255
    pleat[pleat!=50]=0
    darkMo[darkMo!=85]=0
    Mo[Mo!=150]=0
    Pt[Pt!=255]=0
    # mask the original image
    pleatMask = cv2.GaussianBlur(Image, (5,5), 0)
    pleatMask = (pleat/50)*pleatMask
    darkMoMask = cv2.GaussianBlur(Image, (5,5), 0)
    darkMoMask = (darkMo/85)*darkMoMask
    MoMask = cv2.GaussianBlur(Image, (gaussBlur,gaussBlur), 0)
    MoMask = (Mo/150)*MoMask
    PtMask = cv2.GaussianBlur(Image, (gaussBlur,gaussBlur), 0)
    PtMask = (Pt/255)*PtMask
    
    # apply adjusted threshold
    ret,pleat = cv2.threshold(pleatMask, p,255,threshType)
    ret,darkMo = cv2.threshold(darkMoMask, d,255,threshType)
    ret,MoMask = cv2.threshold(MoMask, m,255,threshType)
    ret,PtMask = cv2.threshold(PtMask, pt,255,threshType)
    # put it all together
    threshedImage = np.add(threshedImage, pleat)
    threshedImage = np.add(threshedImage, darkMo)
    threshedImage = np.add(threshedImage, MoMask)
    threshedImage = np.add(threshedImage, PtMask)
    
    if MaskEdges:
    # Impose the masked edge:
        masked,threshMask = maskEdge(ogimage)
        threshedImage = (threshMask/255)*threshedImage
    # Eliminate the black background so only dirt particles/Pt particles are visible
        comb = (blk/255)*threshMask
        threshedImage = threshedImage+cv2.bitwise_not(comb)
    if GetMask:
        return threshedImage,comb
    else:
        return threshedImage

#_____________________________________________________________
##############--MOLYDENUM LOSS APPROXIMATION--################

def isolatePt (image):
# This function filters and thresholds the image in order to estimate the area of exposed Pt.
#  Could use some tweaking. 
	# Apply a bilateral filter to the image to smooth out the texture in Mo and Pt regions:
	# (bilateral filters preserve the edges of each region, which is important in this case)
	poster = makePoster(image)
	# Threshold the image. This is a global threshold. There is probably a better one out there.
	isoPt = regionalThresh(image, poster,p=90,d=120,m=180,pt=180,gaussBlur=3)
	# Make the image into a boolean image
	isoPt = isoPt.astype(np.bool_)
	return isoPt
	
def getFoilArea(img, res=1):
# This is a function based on maskEdge that returns the approximate area of the foil.
# Parameters are:
#	img = nparray of image,
#	res = the resolution of the image, in square microns per pixel
# 	         Default set to 1 micron^2 = 1 pixel
	#Load image
	image = img.copy()
	# Apply a morphological closing to the image to remove small dark areas:
	kernel = np.ones((3,3))
        cimage = cv2.morphologyEx(image, cv2.MORPH_CLOSE, kernel)
	#Threshold and blur image for contour detection
	blur2 = cv2.blur(cimage, (100,100))
	ret,thresh2 = cv2.threshold(blur2,50,255,cv2.THRESH_BINARY)
        thresh2 = cv2.copyMakeBorder(thresh2,5,5,5,5,cv2.BORDER_CONSTANT,value=0)
	#Find the contours of the foil
	can1 = cv2.Canny(thresh2, 50, 200)
	contours,hierarchy = cv2.findContours(can1.copy(),cv2.RETR_TREE,cv2.CHAIN_APPROX_NONE)
	#Isolate and draw the largest contour, which should be the outline of the foil
	#Sort contours by size:
	contours = sorted(contours, key = cv2.contourArea, reverse = True)
	#Largest contour should be the outer edge of foil
	outercnt = contours[0]
	# Foil area is approximately the area of the outer contour multiplied by
	#  the resolution of the image. 
	foilarea = cv2.contourArea(outercnt)*res

	return foilarea
	

def calcPtArea(img, comp=True):
# This function calculates the area of Pt in a binary thresholded image of the foil.
# Calculate WHITE or BLACK area in a binary (boolean) image by counting pixels.
    if str(img.dtype) == 'bool':
        if comp == False:
            inv = np.logical_not(img)
            return np.sum(inv)
        else:
            return np.sum(img)
    else:
        print "Only accept binary images!"

#_____________________________________________________________
##################--MAKING OUTPUT IMAGES--####################

	
#_____________________________________________________________
#####################--MISCELLANEOUS--########################

def show(*images, **kwargs):
# This is a useful function for plotting and displaying any number of images quickly 
# color is the color option ("gray", "jet" etc.) used in plt.imshow
# Rows is the number of rows of subplots of images
# Any number of images can be inserted
        # Get key word arguments:
        color = kwargs.get('color',"gray")
        rows = kwargs.get('rows',1)
	numImgs = len(images)
	# Columns of images is equal to the number of images divided by the specified number 
	# of rows + the modulus of numImgs and rows.
	columns = numImgs/rows + numImgs%rows
	pos = range(len(images))
	for i in pos:
		plt.subplot(rows,columns,i+1),plt.imshow(images[i],color)
	plt.show()
def loadImg (path, flag=cv2.CV_LOAD_IMAGE_GRAYSCALE):
# This is a function for loading images. It basically just solves an issue with 
# cv2.imread and raises an exception if the path is wrong and cv2 returns a NoneType 
# rather than a numpy array.
    image = cv2.imread(path, flag)
    if isinstance(image, type(None)):
    # check if image is an instance of type 'NoneType'
        if not(os.path.exists("path")):
        # Check if the path exists.
            raise Exception("Image file path does not exist.")
            return
        else:
        # This shouldn't come up often, but if the file exists but cv2.imread still
        # spits out a NoneType object, then this should catch the problem.
            raise Exception("Image file path exists, but cv2.imread couldn't load it.")
            return
    else:
        return image