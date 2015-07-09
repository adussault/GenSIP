# This is the main script for running GenSIP on a folder of before and after images
#=========================================================================================
#=====================================###SCRIPT###========================================
#=========================================================================================

'''
# ANALYZING A CLEAN TEST USING THE SAMPLE SET STRING

from GenSIP.cleantests import dirt, mo

sss = "TEST" # sample set string
res = 16 # reslution of image in microns^2 per pixel area
mo.analyzeMoly(sss, res) 
dirt.analyzeDirt(sss, res)

'''
#==============================================================================#
'''
# ANALYZING A PANORAMA IMAGE AND SUBIMAGES
from GenSIP.bigscans import bigfoils

# ANALYZE FROM A PANORAMA IMAGE
panoramaPath = "InputPicts/Foilscans/panos_masks/40360_Q2_panorama.png"
maskPath = "InputPicts/Foilscans/panos_masks/40360_Q2_mask.png"
res = 16
foilname = "40360_2"

bigfoils.analyzePano(panoramaPath, maskPath, 
                     res, foilname, 
                     Quarter="Q2", 
                     MoDirt="Mo", 
                     GenPoster=False, 
                     verbose=False)

# ANALYZE FROM A FOLDER OF SUBIMAGES (requires running analyzePano the first time)
panoramaFolder = "InputPicts/Foilscans/40360_2/sub_imgs_Q2/"
maskFolder = "InputPicts/Foilscans/40360_2/sub_imgs_Q2_mask/"
res = 16
foilname = "40360_2"
bigfoils.analyzeSubImages(panoramaPath, maskPath, 
                          res, foilname, 
                          Quarter="Q2", 
                          MoDirt="Mo", 
                          GenPoster=False, 
                          verbose=False)
'''
#==============================================================================#
'''
# ANALYZING A FOLDER OF IMAGES
from GenSIP import nexus
# Analyzing foils
path = "InputPicts/TEST_Foils"
res = 16
nexus.animorf(path, res, 
                MoDirt='dirt',
                Mask=0, 
                method='cleantests')
                         
nexus.animorf(path, res, 
                MoDirt='dirt',
                Mask=0, 
                method='bigfoils', 
                autoMaskEdges=True)
                         
nexus.animorf(path, res, 
                MoDirt='dirt',
                Mask=0, 
                method='histogram', 
                autoMaskEdges=True)
                         
nexus.animorf(path, res, 
                MoDirt='mo',
                Mask=0,
                method='cleantests')

# ANALYZING SUBIMAGES
path = "InputPicts/TEST_subImages"
masks = "InputPicts/TEST_subImages_masks"

# TEST WITHOUT MASK:
nexus.animorf(path, res, 
                MoDirt='dirt',
                Mask=0, 
                method='histogram', 
                autoMaskEdges=True)
# This will be innacurate b/c of the automatic edge masking feature of cleantests
nexus.animorf(path, res, 
                MoDirt='dirt',
                Mask=0, 
                method='cleantests', 
                autoMaskEdges=True)

# TEST WITH MASK PROVIDED:
nexus.animorf(path, res, 
              MoDirt='dirt',
              Mask=masks, 
              method='histogram', 
              autoMaskEdges=True)
# This will be innacurate b/c of the automatic edge masking feature of cleantests
nexus.animorf(path, res, 
              MoDirt='dirt',
              Mask=masks, 
              method='cleantests', 
              autoMaskEdges=True)
'''
#==============================================================================#
'''
# TEST AN INDIVIDUAL IMAGE

from GenSIP import nexus

path = "InputPicts/TEST_subImages/sub_010_001.tif"
mask = "InputPicts/TEST_subImages_masks/sub_010_001.tif"
res = 16

nexus.animorf(path, res, 
                         MoDirt='dirt',
                         Mask=mask, 
                         method='histogram')
'''