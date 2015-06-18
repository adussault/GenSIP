# This is the main script for running GenSIP on a folder of before and after images
#=========================================================================================
#=====================================###SCRIPT###========================================
#=========================================================================================
import GenSIP_v4.cleantests.moly as mo
import GenSIP_v4.cleantests.dirt as dirt
import GenSIP_v4.functions as fun

sss = "TEST"
mo.analyzeMoly(sss)
dirt.analyzeDirt(sss)