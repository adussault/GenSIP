# This is the main script for running GenSIP on a folder of before and after images
#=========================================================================================
#=====================================###SCRIPT###========================================
#=========================================================================================
import GenSIP.cleantests.moly as mo
import GenSIP.cleantests.dirt as dirt
import GenSIP.functions as fun

sss = "TEST"
mo.analyzeMoly(sss)
dirt.analyzeDirt(sss)