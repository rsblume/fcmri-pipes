#!/usr/bin/python

import os
from os.path import exists
import commands
import sys
import utils
import time
from optparse import OptionParser

class Tcolor:
    HEADER = '\033[0;40;33m';  BLUE = '\033[94m';    ERROR = '\033[0;37;41m'; 
    SUCCESS = '\033[0;42;37m'; WARNING = '\033[93m'; END = '\033[0m'

def print_out(self, text):
	print "\n" + Tcolor.HEADER + text + '\033[0;40;32m \nSubject: ' + self.options.ss  + ", ("+time.ctime()+")" + Tcolor.END

masks= ["WM","VENT"]

class MeanPscExtent:

    def get_opts(self):
        desc="""This is a program to get mean, percent signal change and extent mask for functional time series."""     
        self.usage = "usage: %prog [options]"
        self.parser = OptionParser(description=desc, version="%prog 19.December.2012")
        self.parser.add_option("--subject", dest="ss",
            help="subject identifier-- must match label used for subject directory and HEAD/BRIK files")
        self.parser.add_option("--newnames", dest="nn",
            help="functional scans from individual pathfiles")
        self.parser.add_option("--blur", dest="blur",
            help="amount to ''spatially smooth'' -- specified in individual pathfiles")
        self.parser.add_option("--path", dest="path",
            help="path to preprocessing directory")
        (self.options, args) = self.parser.parse_args()

    def Mean(self):
        os.chdir(self.options.path)
        print "\nCalculating the voxel-wise mean for each functional run for subject "+self.options.ss+" - "+time.ctime()
        for i in self.options.nn.split('\t')[:]:
            idname = self.options.ss+"."+i
#in next line, "orig.BRIK[5..$]" discards first 5 crummy scans...check to see if increasing DS from 2-7 fixes this.
            commands.getoutput("3dTstat -mean -stdev -prefix mean.despike."+idname+" 'reg.despike.shift."+idname+"+orig.BRIK[5..$]'")

    def PerSigChange(self):
        os.chdir(self.options.path)
        print "\nCalculating the percent signal change for each functional run for subject "+self.options.ss+" - "+time.ctime()
        for i in self.options.nn.split('\t')[:]:
            idname = self.options.ss+"."+i
            commands.getoutput("3dcalc -datum float -a reg.despike.shift."+idname+"+orig -b mean.despike."+idname+"+orig[0] -expr 'max(min((a-b)/b*100,100),-100)' -fscale -prefix PSC."+idname)

    def BlurPSC(self):
        os.chdir(self.options.path)
        print "\n''Spatially smoothing'' (um, blurring) each PSC functional run for subject "+self.options.ss+" - "+time.ctime()
        for i in self.options.nn.split('\t')[:]:
            idname = self.options.ss+"."+i
            commands.getoutput("3dmerge -1blur_fwhm "+self.options.blur+" -doall -prefix blur_PSC."+idname+" PSC."+idname+"+orig.")

    def BlurDetRegDesShf(self):
        os.chdir(self.options.path)
        print "\n''Spatially smoothing'' (um, blurring) each detrended, registered, despiked, and time shifted functional run for subject "+self.options.ss+" - "+time.ctime()
        for i in self.options.nn.split('\t')[:]:
            idname = self.options.ss+"."+i
            commands.getoutput("3dmerge -1blur_fwhm "+self.options.blur+" -doall -prefix blur_det.reg.despike.shift."+idname+" det.reg.despike.shift."+idname+"+orig.")

    def MaskExtent (self):
### COMMENTED OUT RIGHT NOW......DO YOU WANT A MASK FOR DECONVOLVE?? (also needs to be fixed if it's a-gonna be used.)
#        self.parser.add_option("--mask", dest="mask",
#            help="appropriate automask dilation from pathfiles") ##  need to fix for each fmri run!! different voxel sizes!
        os.chdir(self.options.path)
        #print os.getcwd()
        for i in self.options.nn.split('\t')[:]:
            idname = self.options.ss+"."+i
#            commands.getoutput("3dcalc -datum float -a reg.despike."+idname+"+orig -expr a -prefix despike."+idname+".float")
            commands.getoutput("3dcalc -datum short -a reg.despike.shift."+idname+"+orig -expr a -prefix reg.despike.shift."+idname+".short")#check input same as 7.mean.py script
            commands.getoutput("3dTstat -datum short -min -prefix mask.extent."+idname+".nostep reg.despike.shift."+idname+".short+orig")
            commands.getoutput("3dcalc -a mask.extent."+idname+".nostep+orig. -b automask."+self.options.ss+"."+self.options.mask+"+orig. -expr '(100+a+b)*step(a+b)' -prefix mask.extent."+idname+".noauto")#check automask dilation
            commands.getoutput("3dAutomask -clfrac 0.1 -prefix mask.extent."+idname+" mask.extent."+idname+".noauto+orig.")
            commands.getoutput("rm -f *noauto* *nostep* reg.despike.shift."+idname+".short")

    def SNR_old (self): 
        os.chdir(self.options.path)
        if exists(self.options.path+"/SNR")==False:
            commands.getoutput("mkdir SNR")      
        ### MAKE directory "SNR" if not existing and put 3 SNR outputs in there: 
        ## one as below, one before despiking (after registering), one mean(despiked, not detrended)/std(despiked,not detrended) [ignore first 5 TRs until DS save the world].
        print "\nMaking SNR directory and populating with signal noise maps by condition for subject "+self.options.ss+"." 
        for i in self.options.nn.split('\t')[:]:
            idname = self.options.ss+"."+i
            commands.getoutput("3dDetrend -polort 3 -prefix det.reg.despike.shift."+idname+" reg.despike.shift."+idname+"+orig.")
            commands.getoutput("3dTstat -mean -stdev -prefix mean.det."+idname+" det.reg.despike.shift."+idname+"+orig.BRIK[5..$]")
            commands.getoutput("3dTstat -mean -stdev -prefix mean.Raw."+idname+" 'reg.shift."+idname+"+orig.BRIK[5..$]'")
            commands.getoutput("3dTstat -mean -stdev -prefix mean.PSC."+idname+" 'PSC."+idname+"+orig.BRIK[5..$]'")
            commands.getoutput("3dTstat -mean -stdev -prefix mean.blur_PSC."+idname+" 'blur_PSC."+idname+"+orig.BRIK[5..$]'")

            commands.getoutput("3dcalc -a mean.Raw."+idname+"+orig.[0] -b mean.Raw."+idname+"+orig.[1] -expr 'a/b' -prefix SNR.Raw."+idname)
            commands.getoutput("3dcalc -a mean.despike."+idname+"+orig.[0] -b mean.despike."+idname+"+orig.[1] -expr 'a/b' -prefix SNR.DespikeOnly."+idname)
            commands.getoutput("3dcalc -a mean.despike."+idname+"+orig.[0] -b mean.det."+idname+"+orig.[1] -expr 'a/b' -prefix SNR.DespikedDetrended_stdev."+idname)
            commands.getoutput("3dcalc -a mean.PSC."+idname+"+orig.[0] -b mean.PSC."+idname+"+orig.[1] -expr 'a/b' -prefix SNR.PSC."+idname)
            commands.getoutput("3dcalc -a mean.blur_PSC."+idname+"+orig.[0] -b mean.blur_PSC."+idname+"+orig.[1] -expr 'a/b' -prefix SNR.blur_PSC."+idname)
        commands.getoutput("mv SNR.* SNR")
            #commands.getoutput("3dcalc -a SNR."+idname+"+orig. -b resampled.anat.automask."+self.options.ss+".2+orig. -expr 'a*step(b)' -prefix SNR.automask."+idname)###AUTOMASK PROBLEM???
        #commands.getoutput("rm -rf det.despike.* stdev.*")

    def SNR (self): #This uses mean vals with discarded first 5 sub-bricks
        os.chdir(self.options.path)
        if exists(self.options.path+"/6.SNR")==False:
            commands.getoutput("mkdir 6.SNR")      
        ### MAKE directory "SNR" if not existing and put 3 SNR outputs in there: 
        ## one as below, one before despiking (after registering), one mean(despiked, not detrended)/std(despiked,not detrended) [ignore first 5 TRs until DS save the world].
	print_out(self, "--- Making SNR directory and populating with signal noise maps by condition")
        for i in self.options.nn.split('\t')[:]:
            idname = self.options.ss+"."+i
            commands.getoutput("3dDetrend -polort 3 -prefix 6.SNR/6.det.reg.despike.shift."+idname+" 3.Register_Despike/3.reg.despike.shift."+idname+"+orig.")
	
	    commands.getoutput("3dTstat -mean -stdev -prefix 6.SNR/6.mean.det.reg.despike.shift."+idname+" 6.SNR/6.det.reg.despike.shift."+idname+"+orig.BRIK[5..$]")
            commands.getoutput("3dTstat -mean -stdev -prefix 6.SNR/6.mean.reg.despike.shift."+idname+" 3.Register_Despike/3.reg.despike.shift."+idname+"+orig.BRIK[5..$]")
            commands.getoutput("3dTstat -mean -stdev -prefix 6.SNR/6.mean.despike.shift."+idname+" 3.Register_Despike/3.despike.shift."+idname+"+orig.BRIK[5..$]")
            commands.getoutput("3dTstat -mean -stdev -prefix 6.SNR/6.mean.ts.shift."+idname+" 3.Register_Despike/3.ts.shift."+idname+"+orig.BRIK[5..$]")
            commands.getoutput("3dTstat -mean -stdev -prefix 6.SNR/6.mean.ts."+idname+" 1.ConvertMR/1.ts."+idname+"+orig.BRIK[5..$]")

            commands.getoutput("3dcalc -a 6.SNR/6.mean.det.reg.despike.shift."+idname+"+orig.[0] -b 6.SNR/6.mean.det.reg.despike.shift."+idname+"+orig.[1] -expr 'a/b' -prefix 6.SNR/6.SNR.det.reg.despike.shift."+idname)
            commands.getoutput("3dcalc -a 6.SNR/6.mean.reg.despike.shift."+idname+"+orig.[0] -b 6.SNR/6.mean.reg.despike.shift."+idname+"+orig.[1] -expr 'a/b' -prefix 6.SNR/6.SNR.reg.despike.shift."+idname)
            commands.getoutput("3dcalc -a 6.SNR/6.mean.despike.shift."+idname+"+orig.[0] -b 6.SNR/6.mean.despike.shift."+idname+"+orig.[1] -expr 'a/b' -prefix 6.SNR/6.SNR.despike.shift."+idname)
            commands.getoutput("3dcalc -a 6.SNR/6.mean.ts.shift."+idname+"+orig.[0] -b 6.SNR/6.mean.ts.shift."+idname+"+orig.[1] -expr 'a/b' -prefix 6.SNR/6.SNR.ts.shift."+idname)
            commands.getoutput("3dcalc -a 6.SNR/6.mean.ts."+idname+"+orig.[0] -b 6.SNR/6.mean.ts."+idname+"+orig.[1] -expr 'a/b' -prefix 6.SNR/6.SNR.ts."+idname)

    def MaskAve(self):
        os.chdir(self.options.path)
        for i in self.options.nn.split('\t')[:]:
            idname = self.options.ss+"."+i
            print "\nGetting WM and ventricle mask values to include in 3dDeconvolve for subject "+self.options.ss
            for mm in masks:
                commands.getoutput("3dmaskave -mask Maskmaker/mask."+mm+"clip.frac."+idname+"+orig.BRIK PSC."+idname+"+orig.BRIK > "+mm+".PSC."+idname+".1D")
                commands.getoutput("awk '{print $1}' "+mm+".PSC."+idname+".1D > col."+mm+".PSC."+idname+".1D")

    def CheckAndClean(self):
        os.chdir(self.options.path)
        for i in self.options.nn.split('\t')[:]:
            idname = self.options.ss+"."+i
#            if exists("mean.despike."+idname+"+orig.BRIK")==False:
#                print "\n\nERROR: This script has NOT successfully completed. Voxelwise mean files have not been created. Error occurred for "+idname+".\n\n"
#                quit()   
#            if exists("PSC."+idname+"+orig.BRIK")==False:
#                print "\n\nERROR: This script has NOT successfully completed. Percent signal change files have not been created. Error occurred for "+idname+".\n\n"
#                quit()
            if exists("6.SNR/6.SNR.det.reg.despike.shift."+idname+"+orig.BRIK")==False:
                print "\n" + Tcolor.ERROR + "ERROR: This script has NOT successfully completed. SNR maps have not been constructed. Error occurred for "+idname+"."+ Tcolor.END + "\n"
                quit()
            if exists("6.SNR/6.SNR.reg.despike.shift."+idname+"+orig.BRIK")==False:
                print "\n" + Tcolor.ERROR + "ERROR: This script has NOT successfully completed. SNR maps have not been constructed. Error occurred for "+idname+"."+ Tcolor.END + "\n"
                quit()
            if exists("6.SNR/6.SNR.despike.shift."+idname+"+orig.BRIK")==False:
                print "\n" + Tcolor.ERROR + "ERROR: This script has NOT successfully completed. SNR maps have not been constructed. Error occurred for "+idname+"."+ Tcolor.END + "\n"
                quit()
            if exists("6.SNR/6.SNR.ts.shift."+idname+"+orig.BRIK")==False:
                print "\n" + Tcolor.ERROR + "ERROR: This script has NOT successfully completed. SNR maps have not been constructed. Error occurred for "+idname+"."+ Tcolor.END + "\n"
                quit()
            if exists("6.SNR/6.SNR.ts."+idname+"+orig.BRIK")==False:
                print "\n" + Tcolor.ERROR + "ERROR: This script has NOT successfully completed. SNR maps have not been constructed. Error occurred for "+idname+"."+ Tcolor.END + "\n"
                quit()
#            for mm in masks:
#                if exists("col."+mm+".PSC."+idname+".1D")== False:
#                    print "\n\nERROR: This script has NOT successfully completed. 1D files for WM and ventricle mask values have NOT been generated. Check subject and file names and data location, etc. Error occurred for "+mm+" mask for run "+idname+".\n\n"
#                    quit()                 
        print "\n" + Tcolor.SUCCESS + "Congratulations! This script (6.PSCmeanSNR--really just calculates 5 SNRs) appears to have run successfully! Check that all of the output files have values. Check the SNR maps in the "+self.options.path+"6.SNR directory!" + Tcolor.END + "\n"
      

mpe=MeanPscExtent()
mpe.get_opts()
#mpe.Mean()-----------------don't need this anymore, i do it in SNR
#mpe.PerSigChange()---------never understood the point of this (well, supposedly a normalizing technique...but deconvolve step takes care of potential image scaling difference between subjects
#mpe.BlurPSC()--------------see above...why blur something not useful?
#mpe.BlurDetRegDesShf()-----as of 2014-10-09, blurring occurs in the next step, deconvolve, but potentially should be here.
########mpe.MaskExtent()  DO YOU WANT A MASK FOR DECONVOLVE??---yes, but this needs some work
mpe.SNR()
##mpe.MaskAve()----------needs some work
mpe.CheckAndClean()
