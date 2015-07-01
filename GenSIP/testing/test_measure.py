"""
Performs tests on calcDirt and calcExposedPt using the small test images stored 
in 
"""

import numpy as np
import GenSIP.functions as fun
import GenSIP.measure as meas
import os
import mahotas as mh
import unittest
import nose


"""
_________________________________
CREATING TESTS USING SMALL IMAGES\______________________________________________

"""
class Test_Measure_With_Small_images (unittest.TestCase):
    
    def setUp(self):
        self.DIRNAME = os.path.split(__file__)[0]
        self.BWFolder = os.path.join(self.DIRNAME,'small_imgs_BW')
        
        if not os.path.exists(self.BWFolder):
            raise Exception(str(os.getcwd()) + " Folder does not exist: " + self.BWFolder)
        self.BW_names = [m for m in os.listdir(self.BWFolder) if m.endswith('.png')]
        
        self.BWImgs = {} # Initiate dictionary of Black & White small images
        self.Exp_Res_Sm_Imgs = {} # Initiate dictionary of expected results
        
        for name in self.BW_names:
            
            self.BWImgs[name] = fun.loadImg(os.path.join(self.BWFolder,name))
            self.BWImgs[name][self.BWImgs[name]!=self.BWImgs[name].max()]=0
            self.Exp_Res_Sm_Imgs[name] = {}
            
            self.Exp_Res_Sm_Imgs[name]['whitePxCount'] = \
            self.BWImgs[name][self.BWImgs[name]==self.BWImgs[name].max()].size
            
            expectedLabeled_3x3BC,expectedNum = mh.label(self.BWImgs[name],Bc=np.ones((3,3)))
            expectedSizes = mh.labeled.labeled_size(expectedLabeled_3x3BC)
            expectedSizes = np.sort(expectedSizes)[::-1]
            expectedSizes = expectedSizes[1:]
            
            self.Exp_Res_Sm_Imgs[name]['expectedLabeled_3x3BC'] = expectedLabeled_3x3BC
            self.Exp_Res_Sm_Imgs[name]['expectedNum'] = expectedNum
            self.Exp_Res_Sm_Imgs[name]['expectedSizes'] = expectedSizes
    """
    Tests for the basic calcDirt and calcExposedPt in GenSIP.measure:
    """
    def test_BWFolder_exists(self):
        nose.tools.assert_true(os.path.exists(self.BWFolder), msg=str(os.listdir('.')))

    def test_BWImgs_is_a_dictionary(self):
        nose.tools.assert_equal(type(self.BWImgs),dict,msg=str(self.BWImgs))
        
    def test_BWImgs_is_a_dictionary_of_images(self):
        for img in self.BWImgs:
            nose.tools.assert_equal(type(self.BWImgs[img]),
                                    np.ndarray,
                                    msg=str(type(self.BWImgs[img])))
            
    def test_calcDirt_Num_Accuracy_on_Small_Imgs_No_Minimum_Particle_Area(self):
        """Test calcDirt on the small images: Accuracy of dirt particle number"""
        Dirt_areas = {}
        Dirt_nums = {}
        for img in self.BWImgs:
            Dirt_areas[img],Dirt_nums[img] = meas.calcDirt(self.BWImgs[img],1)
        for img in Dirt_nums:
            nose.tools.assert_equal(Dirt_nums[img],self.Exp_Res_Sm_Imgs\
            [img]['expectedNum'], 
            msg="""
                calcDirt does not get the right area for small image {0}: \n
                Expected Value: {1} \n
                calcDirt Value: {2}
                """.format(
                        img,
                        self.Exp_Res_Sm_Imgs[img]['expectedNum'],
                        Dirt_nums[img]))
            print img
        
    def test_calcDirt_Area_Accuracy_on_Small_Imgs_No_Minimum_Particle_Area(self):
        """Test calcDirt on the small images: Accuracy of dirt area"""
        Dirt_areas = {}
        Dirt_nums = {}
        for img in self.BWImgs:
            Dirt_areas[img],Dirt_nums[img] = meas.calcDirt(self.BWImgs[img],1)
        for img in Dirt_areas:
            nose.tools.assert_equal(Dirt_areas[img], 
            self.Exp_Res_Sm_Imgs[img]['whitePxCount'], 
            msg="""
                calcDirt does not get the right area for small image {0}: \n
                Expected Value: {1} \n
                calcDirt Value: {2}
                """.format(
                        img,
                        self.Exp_Res_Sm_Imgs[img]['whitePxCount'], 
                        Dirt_areas[img]))
            print img
    
    def test_calcDirt_Sizes_Accuracy_on_Small_Imgs_No_Minimum_Particle_Area(self):
        Dirt_sizes = {}
        for img in self.BWImgs:
            Dirt_sizes[img] = meas.calcDirt(self.BWImgs[img],1,returnSizes=True)[2]
        for img in Dirt_sizes:
            nose.tools.assert_equal(
                Dirt_sizes[img].all(), 
                self.Exp_Res_Sm_Imgs[img]['expectedSizes'].all(), 
                msg="""
                    calcDirt does not give the correct sizes for small image {0}: \n
                    Expected Values: {1} \n
                    calcDirt Values: {2}
                    """.format(
                            img,
                            self.Exp_Res_Sm_Imgs[img]['expectedSizes'], 
                            Dirt_sizes[img]))
            print img
            
    def test_calcDirt_Sizes_Accuracy_on_Small_Imgs_10micron_Minimum_Particle_Area(self):
        Dirt_sizes = {}
        for img in self.BWImgs:
            Dirt_sizes[img] = meas.calcDirt(self.BWImgs[img],
                                            1,
                                            minPartArea=10,
                                            returnSizes=True)[2] #sizes is at index 2
        for img in Dirt_sizes:
            Expcted_sizes = self.Exp_Res_Sm_Imgs[img]['expectedSizes'][
                                self.Exp_Res_Sm_Imgs[img]['expectedSizes']>=10]
            nose.tools.assert_equal(
                Dirt_sizes[img].all(), Expcted_sizes.all(), 
                msg="""
                    calcDirt does not give the correct sizes for small image {0}: \n
                    Expected Value: {1} \n
                    calcDirt Value: {2}
                    """.format(
                            img,
                            self.Exp_Res_Sm_Imgs[img]['expectedSizes'], 
                            Dirt_sizes[img]))
            print img
            
    def test_calcDirt_Sizes_on_Small_Imgs_Minimum_Area_Larger_than_Img(self):
        Dirt_sizes = {}
        for img in self.BWImgs:
            minArea = self.BWImgs[img].size+1
            Dirt_sizes[img] = meas.calcDirt(self.BWImgs[img],
                                            1,
                                            minPartArea=minArea,
                                            returnSizes=True)[2] #sizes is at index 2
            Expcted_sizes = np.array([])
            nose.tools.assert_equal(
                Dirt_sizes[img].all(), 
                Expcted_sizes.all(), 
                msg="""
                    calcDirt does not give the correct sizes for small image {0}: \n
                    Expected Value: {1} \n
                    calcDirt Value: {2}
                    """.format(
                            img,
                            self.Exp_Res_Sm_Imgs[img]['expectedSizes'], 
                            Dirt_sizes[img]))
            print img
            
    def test_calcDirt_Labelled_Img_on_Small_Imgs_Bc_3x3_Ones(self):
        """Test to make sure the labelled image comes out correctly, 
           BoundConds = 3x3 array of ones (8 nearest neighbors)"""
        Dirt_labelled_imgs = {}
        for img in self.BWImgs:
            Dirt_labelled_imgs[img] = meas.calcDirt(self.BWImgs[img],
                                                    1,
                                                    returnSizes=False,
                                                    BoundConds=np.ones((3,3)),
                                                    returnLabelled=True)[2]
        for img in Dirt_labelled_imgs:
            nose.tools.assert_equal(
                Dirt_labelled_imgs[img].all(), 
                self.Exp_Res_Sm_Imgs[img]['expectedLabeled_3x3BC'].all(), 
                msg="""
                    calcDirt does not give the correct labelled image for small image {0}: \n
                    Expected Image: \n 
                    {1} \n
                    calcDirt Image: \n 
                    {2}
                    """.format(
                            img,
                            self.Exp_Res_Sm_Imgs[img]['expectedLabeled_3x3BC'], 
                            Dirt_labelled_imgs[img]))
            print img
        
    def test_calcExposedPt_Area_Accuracy_Small_Imgs(self):
        Pt_areas = {}
        for img in self.BWImgs:
            Pt_areas[img] = meas.calcExposedPt(self.BWImgs[img], 1)
        for img in Pt_areas:
            print img + ' passed.'
            print round(Pt_areas[img],8)
            print int(Pt_areas[img]*1000000)
            print self.Exp_Res_Sm_Imgs[img]['whitePxCount']
            nose.tools.assert_almost_equal(Pt_areas[img]*1000000, 
            self.Exp_Res_Sm_Imgs[img]['whitePxCount'], 
            delta=.1,
            msg="""
                calcExposedPt does not get the right area for small image {0}: \n
                Expected Value: {1} \n
                calcExposedPt Value: {2}
                """.format(img, 
                    self.Exp_Res_Sm_Imgs[img]['whitePxCount'], 
                    Pt_areas[img]*1000000))
        
    
if __name__=='__main__':
    import sys
    print os.getcwd()

    module_name = sys.modules[__name__].__file__
    nose.run(argv=[sys.argv[0],module_name,'-v'])
    nose.run(defaultTest=__name__)
