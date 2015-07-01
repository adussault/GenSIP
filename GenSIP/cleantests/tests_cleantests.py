import unittest
import nose
import os
import GenSIP.functions as fun


class Test_Comp_Funcs(unittest.TestCase):
    def setUp(self):
        self.DIRNAME = os.path.split(__file__)[0]
        self.BefPict = fun.loadImg(os.path.join(self.DIRNAME,"40393,0210 before clean.tif"))
        self.AftPict = fun.loadImg(os.path.join(self.DIRNAME,"40393,0210 after clean.tif"))
    def test_