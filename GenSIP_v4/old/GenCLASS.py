import GenSIP_v4.functions as funcs
import cv2
import numpy as np

class ContourRegistry(object):
    def __init__(self,Kimage):
        self.Image = Kimage
        self.Contours,\
        self.Hierarchy,\
        self.Canny = self.getRegions(Kimage)
        self.Regions = self.RegConts(Kimage)
    def RegConts(self,image):
        numConts = len(self.Contours)
        regions = []
        for i in range(numConts):
            regions.append(Region(image,self.Contours,self.Hierarchy,i))
        return regions
    def getRegions(self, regImg):
        #Find the contours of the foil
        can = cv2.Canny(regImg, 50, 200)
        conts,hier = cv2.findContours(can.copy(),cv2.RETR_TREE,cv2.CHAIN_APPROX_NONE)
        return conts,hier,can
    # Not a stable method:
    '''
    def fixAreas(self): 
        for r in self.Regions:
            print "Got to Region %s" %str(r.Index)
            perc = int(r.Area/r.ContourArea*100)
            if r.Area == 0:
                print "Region %s has no area." %str(r.Index)
            # if the area of the region is less than 5 percent of the area
            # of its outline and the area of the region is < 100 pixels,
            # then remove it from the hieararchy and reincorporate its mask
            # into another region (all done in that region's removeHier method
            if perc < 10 and r.Area <=100:
                #print "Region %s area: %s" %(r.Index,r.Area)
                r.removeHier(self)
                print "Region %s removed from hierarchy." %str(r.Index)
            elif r.ContourArea < 20 or r.Area < 20:
               # print "Region %s area: %s" %(r.Index,r.Area)
                r.removeHier(self)
                print "Region %s removed from hierarchy." %str(r.Index)
        print "Resetting Hierarchy..."
        self.resetHier()
    '''
    def fixAreas(self):
        fixes = []
        for r in self.Regions:
            #print "Got to Region %s" %str(r.Index)
            if r.FirstChild != -1 and r.Area < self.Regions[r.FirstChild].Area:
                print "Region %s has a smaller area than child %s" \
                %(r.Index,r.FirstChild)
                print "Area: %s" % r.Area
                r.removeHier()
                print "Region %s removed" % r.Index
                fixes.append((r.Index,r.FirstChild))
        return fixes
    def finalThresh(self,image):
        image_c = image.copy()
        canvas = np.zeros(image.shape,np.uint8)
        for r in self.Regions:
            mask = (r.Mask/255)*image_c
            if r.GrayLevel <= 10:
                ret,thresh = cv2.threshold(mask, 0,255,cv2.THRESH_BINARY)
            elif 40<=r.GrayLevel<=60:
                ret,thresh = cv2.threshold(mask, 0,255,cv2.THRESH_BINARY)
            elif 75<=r.GrayLevel<=100:
                ret,thresh = cv2.threshold(mask, 0,255,cv2.THRESH_BINARY)
            elif 140<=r.GrayLevel<=160:
                ret,thresh = cv2.threshold(mask, 0,255,cv2.THRESH_BINARY)
            elif r.GrayLevel>220:
                ret,thresh = cv2.threshold(mask, 0,255,cv2.THRESH_BINARY)
            else: 
                print "Region %s falls out of range." % r.Index
            canvas = np.add(canvas,thresh)
        return canvas

    def makeMask(self,regIndex):
    # This method creates an image of the original image masked by this region. 
        maskedImage = self.Image.copy()
        mask = (self.Regions[regIndex].Mask==0)
        maskedImage[mask] = 0
        return maskedImage
    def resetHier(self):
    # This resets the hierarchy by looking at the Nodes of each region.
        for r in self.Regions:
            self.Hierarchy[r][0]=r.Next
            self.Hierarchy[r][1]=r.Prev
            self.Hierarchy[r][2]=r.FirstChild
            self.Hierarchy[r][3]=r.Parent
            # also reset the nodes to make sure they correspond to the 
            # rest of the attributes
            r.resetNode()
        print "Hierarchy reset."
    
class Region(object):
    def __init__(self,image,conts,hier,index):
        # Image is a Kuwahara and posterized image with the same resolution as
        #  the original foil image.
        #self.Image = image #May want to make a foil image class in the future to store the resolution, foilname, size, Pt area, Dirtcount, etc. 
        self.ImageShape = image.shape
        # Index is the index of the region in the contour list
        self.Index = index
        self.Contour = conts[index]
        # Set up Hierarchy relationships:
        self.Node = hier[0,index] #4-length array defining place in hierarchy tree
        self.Next = self.Node[0]
        self.Prev = self.Node[1]
        self.FirstChild = self.Node[2]
        self.Parent = self.Node[3]
        # Create mask
        self.Outline = self.drawOutline(conts,hier)
        self.Mask = self.maskCont(self.ImageShape,conts,hier,index)
        # Region properties
        self.ContourArea=cv2.contourArea(self.Contour)# May want to mult by resolution in future
        self.OutlineArea=np.sum(self.Outline)/255
        self.Area = np.sum(self.Mask)/255
        # Mean gray level
        #cv2 mean funtion takes an image and a mask as arguments and 
        # returns an array of the mean pixel values as floats. We only need the first value
        self.GrayLevel = int(cv2.mean(image,self.Mask)[0])
        
    def drawOutline(self,conts,hier):
        canvas = np.zeros(self.ImageShape,np.uint8)
        cv2.drawContours(canvas,conts,self.Index,255,-1)
        return canvas
        
    def maskCont(self,shape,conts,hier,index):
    # Creates a mask for a single contour, cutting out any contours that might lie
    #     within that contour. Arguments are the image shape tuple, the contours in
    #     image, the hierarchy tree, and the index of the mask to be masked.
        canvas = self.Outline.copy()
        children = np.zeros(shape,np.uint8)
        # Check if there is a first child:
        if self.FirstChild!=-1:
            children=self.drawAndSiblings(children,conts,hier,self.FirstChild)
        # Subtract out the children
        canvas = canvas-children
        return canvas
        
    def drawAndSiblings(self,canvas, conts, hier, index):
    # This method is a tool for drawing the region and its siblings. It is mainly used
    # in the maskCont function in order to ensure that all child areas are excluded
    # from the mask of the parent region.
        nxt = hier[0,index,0]
        cv2.drawContours(canvas,conts,index,255,-1)
        if nxt!=-1:
            canvas=self.drawAndSiblings(canvas,conts,hier,nxt)
        return canvas
        
    
    
    def removeHier(self):
        # takes a ContourRegistry object (Regs) and removes this contour from the
        # hierarchy, then tells Regs to reconstruct its hierarchy and deletes self.
        # Remove from Siblings
        if self.Next != -1:
            Regs.Regions[self.Next].Prev = self.Prev
        if self.Prev != -1:
            Regs.Regions[self.Prev].Next=self.Next
        # if there is a first child, set the parent of the first child to the 
        # parent of self, then set the previous of first child to prev of self, 
        # then iterate through the siblings of firstchild and set the parents to
        # the parent of self. When there are no more next children in line, set 
        # the last child's next to the next of self, and previous of next to last child
        # and then set the parent of all 
        if self.FirstChild != -1:
            Regs.Regions[self.FirstChild].Parent = self.Parent
            Regs.Regions[self.FirstChild].Prev = self.Prev
            Regs.Regions[self.Prev].Next = self.FirstChild
            i=Regs.Regions[self.FirstChild].Next
            j=Regs.Regions[self.FirstChild].Index
            print "%s while loop started." %self.Index
            print "Firstchild %s next value: %s" %(self.FirstChild, i)
            while i != -1:
                Regs.Regions[i].Parent=self.Parent
                j = Regs.Regions[i].Index
                i=Regs.Regions[i].Next
            print "%s while loop finished." %self.Index
            Regs.Regions[j].Next = self.Next
            if self.Next != -1:
                Regs.Regions[self.Next].Prev = j
            
        # if there is a parent and self is the first child of the parent, then
        #  set the first child of the parent to the first child of self.
        if self.Parent != -1 and Regs.Regions[self.Parent].FirstChild==self.Index:
            Regs.Regions[self.Parent].FirstChild = self.FirstChild
             
               
        # Incorporate pixels into another mask. Incorporate it into the darker mask
        # for optimal thresholding:
        if self.Parent != -1 and self.FirstChild != -1:
        # If there is both a parent and a child, check to see which has a lower mean GrayLevel
        # and then incorporate the pixels from this region into that region's mask.
            if Regs.Regions[self.Parent].GrayLevel <= Regs.Regions[self.FirstChild].GrayLevel:
                Regs.Regions[self.Parent].Mask = np.add(self.Mask,Regs.Regions[self.Parent].Mask)
                print "Region %s pixels incorporated to parent region %s." %(str(self.Index),str(self.Parent))
            else:
                Regs.Regions[self.FirstChild].Mask = np.add(self.Mask,Regs.Regions[self.FirstChild].Mask)
                print "Region %s pixels incorporated to child region %s." %(str(self.Index),str(self.FirstChild))
        elif self.Parent != -1:
        # If only a parent exists, incorporate the mask into the parent's mask
            Regs.Regions[self.Parent].Mask = np.add(self.Mask,Regs.Regions[self.Parent].Mask)
            print "Region %s pixels incorporated to parent region %s.B" %(str(self.Index),str(self.Parent))
        elif self.FirstChild != -1:
        # If only a first child exists, incorporate the mask into the first child's mask
            Regs.Regions[self.FirstChild].Mask = np.add(self.Mask,Regs.Regions[self.FirstChild].Mask)
            print "Region %s pixels incorporated to child region %s. B" %(str(self.Index),str(self.FirstChild))
        else:
            print "Region %s has a tiny area and no parents/children. " %str(self.Index)
        
        # set all neighbors to -1 to remove self from hierarchy
        self.Next = -1
        self.Prev = -1
        self.FirstChild = -1
        self.Parent = -1
        self.resetNode()
        
        #Tell the Contour Registry to reorganize the hierarchy.
        # Should only be done after a series of removals to save time, rather
        # than resetting after each region is removed.
        #Regs.resetHier()
        
    def resetNode(self):
        # resets the node attribute according to the Next, Prev, FirstChild, and 
        # Parent attributes.
        self.Node[0]=self.Prev
        self.Node[1]=self.Next
        self.Node[2]=self.FirstChild
        self.Node[3]=self.Parent

b23 = funcs.loadImg("InputPicts/Before/Before_TEST/40393,0230 before clean.tif",cv2.CV_LOAD_IMAGE_GRAYSCALE)
k23 = funcs.loadImg("InputPicts/Before/Before_TEST/0231 Resize Kuwahara.tif",cv2.CV_LOAD_IMAGE_GRAYSCALE)
mk23 = funcs.maskEdge(k23,thickness=30)
mb23 = funcs.maskEdge(b23,thickness=30)

p23=funcs.posterfy(mk23)
global Regs
Regs=ContourRegistry(p23)
regs=Regs.Regions