import os
#from os.path import exists
#import commands
#import sys
#import utils

class SkullZero:

    def build_skullzero_dirs(self, subjectdir):
	""" This method builds the directory for skullstripping and zeropadding"""
	os.chdir(subjectdir + "/anatomy/")
        command= "mkdir 2.Skullstrip_Zeropad"
	os.system(command)
    
    def AlignAnatCenter(self, subjectdir, anatbasename):
	""" This method copies anatomy puts into folder and then runs 3drefit to center the anatomical image at the origin"""
	command= """cp %s %s/anatomy/2.Skullstrip_Zeropad""" %(anatbasename, subjectdir)
	os.system(command)
	os.chdir(subjectdir + "/anatomy/2.Skullstrip_Zeropad")
	command = """3drefit -xorigin cen -yorigin cen -zorigin cen %s/anatomy/%s """ %(subjectdir, anatbasename)
	os.system(command)


    def SkullStrip(self,subjectdir,anatbasename):
        """ This method skull strips ananatomical image"""
	command = """3dSkullStrip -input %s/anatomy/%s -prefix %s/anatomy/2.Skullstrip_Zeropad/2.volume.stripped.%s""" %(subjectdir,anatbasename,subjectdir, anatbasename)
        os.system(command)


    def ZeroPad(self,padvals, subjectdir,anatbasename):
        """ This method zero-pads the ananatomical image"""
	command = """3dZeropad %s -overwrite -prefix %s/anatomy/2.Skullstrip_Zeropad/2.volume.padded.stripped.%s %s/anatomy/2.Skullstrip_Zeropad/2.volume.stripped.%s""" %(padvals, subjectdir,anatbasename, subjectdir, anatbasename)
        os.system(command)
	#I am no zero padding the result of 1 because I do not have athis step



subjectdir = "/cnari/sports_2014/Data/Subjects/s01_sports_2014"
anatbasename = "s01_sports_2014_WIP_RestSlab_T1W_3D_TFE_SENSE_3_1.nii"
skullstripname =  "2.volume.stripped"
padvals= "-I 20 -S 20 -L 10 -R 10 -A 30 -P 30"
