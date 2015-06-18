################################################################################
###     BAR is a temporary file for writing and testing new functions        ###
################################################################################

execfile("GenSIP/GenSIP.py")
execfile("GenSIP/Kuwahara.py")
execfile("GenSIP/foo.py")

    
    
import cv2 
from scipy import misc
b228 = loadImg("InputPicts/Before/Before_TEST/40393,0228 before clean.tif",cv2.CV_LOAD_IMAGE_GRAYSCALE)
b213 = loadImg("InputPicts/Before/Before_TEST/40393,0213 before clean.tif",cv2.CV_LOAD_IMAGE_GRAYSCALE)
b230 = loadImg("InputPicts/Before/Before_TEST/40393,0230 before clean.tif",cv2.CV_LOAD_IMAGE_GRAYSCALE)
b231 = loadImg("InputPicts/Before/Before_TEST/40393,0231 before clean.tif",cv2.CV_LOAD_IMAGE_GRAYSCALE)
b232 = loadImg("InputPicts/Before/Before_TEST/40393,0232 before clean.tif",cv2.CV_LOAD_IMAGE_GRAYSCALE)

print "Reducing image to 1/10th the size..."
#resize image to 10% of original
rsz228 = misc.imresize(b228,.1,interp='bicubic')
rsz213 = misc.imresize(b213,.1,interp='bicubic')
rsz230 = misc.imresize(b230,.1,interp='bicubic')
rsz231 = misc.imresize(b231,.1,interp='bicubic')
rsz232 = misc.imresize(b232,.1,interp='bicubic')

print "Applying first Gaussian blur..."
# apply Gaussian filter
gaur228 = cv2.GaussianBlur(rsz228, (3,3),0)
gaur213 = cv2.GaussianBlur(rsz213, (3,3),0)
gaur230 = cv2.GaussianBlur(rsz230, (3,3),0)
gaur231 = cv2.GaussianBlur(rsz231, (3,3),0)
gaur232 = cv2.GaussianBlur(rsz232, (3,3),0)

print "Applying Kuwahara filter..."
# apply Kuwahara Filter:
kuwgr228 = Kuwahara(gaur228,13)
kuwgr213 = Kuwahara(gaur213,13)
kuwgr230 = Kuwahara(gaur230,13)
kuwgr231 = Kuwahara(gaur231,13)
kuwgr232 = Kuwahara(gaur232,13)

print "Resizing to original size..."
#resize images to original size
rkgr228 = misc.imresize(kuwgr228,(b228.shape),interp='bilinear')
rkgr213 = misc.imresize(kuwgr213,(b213.shape),interp='bilinear')
rkgr230 = misc.imresize(kuwgr230,(b230.shape),interp='bilinear')
rkgr231 = misc.imresize(kuwgr231,(b231.shape),interp='bilinear')
rkgr232 = misc.imresize(kuwgr232,(b232.shape),interp='bilinear')

print "Applying second Gaussian blur..."
# Gauss blur one more time
grkgr228 = cv2.GaussianBlur(rkgr228, (11,11),0)
grkgr213 = cv2.GaussianBlur(rkgr213, (11,11),0)
grkgr230 = cv2.GaussianBlur(rkgr230, (11,11),0)
grkgr231 = cv2.GaussianBlur(rkgr231, (11,11),0)
grkgr232 = cv2.GaussianBlur(rkgr232, (11,11),0)

'''
# Kuwahara filter once more
kgrkgr230 = Kuwahara(grkgr230,13)
kgrkgr231 = Kuwahara(grkgr231,13)
kgrkgr232 = Kuwahara(grkgr232,13)
# posterfy the image
b230,b231,b232
'''
print "Posterizing..."
pp228 = posterfy(rkgr228)
pp213 = posterfy(rkgr213)
pp230 = posterfy(rkgr230)
pp231 = posterfy(rkgr231)
pp232 = posterfy(rkgr232)
print "Done, Displaying..."
show(grkgr230,grkgr231,grkgr232,pp230,pp231,pp232,rows=2)
print "Display closed."

def proc4Dirt(image):
    rsz = misc.imresize(image,.1,interp='bicubic')
    gr = cv2.GaussianBlur(rsz, (3,3),0)
    kgr = Kuwahara(gr,13)
    rkgr = misc.imresize(kgr,(image.shape),interp='bicubic')
    grkgr = cv2.GaussianBlur(rkgr, (11,11),0)
    prkgr = posterfy(grkgr)
    return prkgr