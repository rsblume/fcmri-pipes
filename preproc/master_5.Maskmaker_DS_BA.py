#!/usr/bin/python
#NOTE: CHECK BLUR in line 78 and in "+ss+".around.vent.blur!!!!

import os
from os.path import exists
import commands
import sys
import utils
import time
from optparse import OptionParser


class MaskMaker:

    def get_opts(self):
        desc="""This is a program to make white matter and ventricle masks to then regress out some sources of nuisance noise.  It uses AFNI functions, as well as FSL.   Note: Make sure to include '+orig' for all input AFNI data types specified here."""      
        self.usage = "usage: %prog [options]"
        self.parser = OptionParser(description=desc, version="%prog 28.April.2010")
        self.parser.add_option("--identity", dest="id",
            help="subject identifier-- must match label used for subject directory and HEAD/BRIK files")
        self.parser.add_option("--volume", dest="vol",
            help="this is your volume data- should be obliqued and aligned T1 anatomical")
        self.parser.add_option("--newnames", dest="nn",
            help="suffix for files for individual runs")
        self.parser.add_option("--location", dest="loc",
            help="Specify the working directory. this option is for the benefit of swift functionality.")
        self.parser.add_option("--path", dest="path",
            help="Path to preprocessed data")
        self.parser.add_option("--mask", dest="mask",
            help="Automask to be used -- specified in pathfiles")               
        (self.options, args) = self.parser.parse_args()
        if len(args) !=0:
            self.parser.error("your arguments == NO BUENO! maybe you put in a flag without giving the argument? maybe you're just messing with me??")
              
    def TLRC_algn(self):
        print """Zeropad and get a tal align coordinate vector"""
        #print "vol= "+self.options.vol
        #print "id= "+self.options.id
        #print "nn= "+self.options.nn
        #print "loc= "+self.options.loc
        #print "mask= "+self.options.mask

        ###??? rm zeropad???
        
        print "<<<<<<<<<<<<<<<<<<<<<< 3dZeropad >>>>>>>>>>>>>>>>>>>>>>>>>>> "+time.ctime()
        commands.getoutput("3dZeropad -I 10 -S 10 -A 10 -P 10 -L 10 -R 10 -prefix "+self.options.loc+"/volume.padded."+self.options.id+" "+self.options.vol) ## create extra padding in the volume space### 
        print "<<<<<<<<<<<<<<<<<<<<<< 3dresample >>>>>>>>>>>>>>>>>>>>>>>>>>> "+time.ctime()
        commands.getoutput("3dresample -inset /soft/afni/TT_N27+tlrc -master "+self.options.loc+"/volume.padded."+self.options.id+"+orig -prefix "+self.options.loc+"/template.resampled."+self.options.id+"+orig") ## resample to tal
        print "<<<<<<<<<<<<<<<<<<<<<< 3dWarpDrive >>>>>>>>>>>>>>>>>>>>>>>>>>> "+time.ctime()
        commands.getoutput("3dWarpDrive -affine_general -cubic -input "+self.options.loc+"/volume.padded."+self.options.id+"+orig -prefix "+self.options.loc+"/tal.volume.padded."+self.options.id+" -base "+self.options.loc+"/template.resampled."+self.options.id+"+orig -1Dmatrix_save "+self.options.loc+"/tal2volume."+self.options.id+".1D") ## using this to get the tal transform matrix 'tal2volume.*.1D'

    def CalcMaskedVol(self):
        print """Get masked data from the volume for ventricles and white matter."""
        auto_mask = self.options.mask
        frac_mask = self.options.loc+"/"+self.options.id+".automask_frac+orig"
        print "<<<<<<<<<<<<<<<<<<<<<< 3dfractionize >>>>>>>>>>>>>>>>>>>>>>>>>>> "+time.ctime()
        commands.getoutput("3dfractionize -template "+self.options.vol+" -input "+auto_mask+" -preserve -clip 0.2 -prefix "+frac_mask)
        print "Sourcing FSLDIR:"
        commands.getoutput("source /soft/fsl/etc/fslconf/fsl.sh") ## source FSL
        print "<<<<<<<< 3dcalc to get masked volume from the fractionized and boxed volume >>>>>>>>>>>>> "+time.ctime()
        commands.getoutput("3dcalc -a "+frac_mask+" -b "+self.options.vol+" -expr 'step(a)*b' -prefix "+self.options.loc+"/volume.masked."+self.options.id+"+orig") ## mask sample volume
        print "<<<<<<<< 3dAFNItoNIFTI to get nii volume >>>>>>>>>>>>> "+time.ctime()
        commands.getoutput("3dAFNItoNIFTI -prefix "+self.options.loc+"/volume.masked."+self.options.id+".nii "+self.options.loc+"/volume.masked."+self.options.id+"+orig") ## convert to nii
        print "<<<<<<<< running 'fast' to get segmentation (not actually fast-- kinda slow) >>>>>>>>>>>>> "+time.ctime()
        os.system(". ${FSLDIR}/etc/fslconf/fsl.sh")
        os.system("fast -o "+self.options.loc+"/volume."+self.options.id+" "+self.options.loc+"/volume.masked."+self.options.id+".nii") ## segments volume
        print "<<<<<<<< 3dcalc to get mask around ventricles from nii segmentation >>>>>>>>>>>>> "+time.ctime()
        commands.getoutput("3dcalc -a "+self.options.loc+"/volume."+self.options.id+"_seg.nii.gz -b "+frac_mask+" -expr '100*step(b)*iszero(amongst(a,1,0))' -prefix "+self.options.loc+"/"+self.options.id+".around.vent.preblur")## mask around the ventricle
        print "<<<<<<<< 3dmerge to blur aroundVent >>>>>>>>>>>>> "+time.ctime()
        ## the blur is not 'one size fits all', i.e., may have to change 1filter_nzmean value so it doesn't blur over the ventricle completely.  start with '-1filter_nzmean 3' per Mike Andric but using 2 here.... seems to work
        commands.getoutput("3dmerge -1filter_nzmean 2 -prefix "+self.options.loc+"/"+self.options.id+".around.vent.blur "+self.options.loc+"/"+self.options.id+".around.vent.preblur+orig")## blur around the ventricle

    def CalcMaskedVolVent(self):
        print "<<<<<<<< Generating the Ventricles seed 'Vseed' >>>>>>>>>>>>> "+time.ctime()
        Vseed = "-8 13 19\n8 13 19\n" ## creating seed in areas of ventricle - tal coords (may not be best coords for all- something to check if output is not existent or ventricles incluse area outside skull)
        file = open(self.options.loc+"/vent.seed.tal.1D",'w')
        file.write(Vseed)
        file.close()
        auto_mask = self.options.mask
        frac_mask = self.options.loc+"/"+self.options.id+".automask_frac+orig"
        print "<<<<<<<< Warping from talairach using 'Vecwarp'  >>>>>>>>>>>>> "+time.ctime()
        os.system("Vecwarp -matvec "+self.options.loc+"/tal2volume."+self.options.id+".1D -forward -input "+self.options.loc+"/vent.seed.tal.1D -output "+self.options.loc+"/vent.seed."+self.options.id+".1D")
        print "<<<<<<<< 3dUndump to get ventricles seed >>>>>>>>>>>>> "+time.ctime()
        commands.getoutput("3dUndump -xyz -orient RAI -prefix "+self.options.loc+"/vent.seed."+self.options.id+" -master "+frac_mask+" -srad 8 "+self.options.loc+"/vent.seed."+self.options.id+".1D")
        print "<<<<<<<< 3dcalc to get inverted vent volume seed  >>>>>>>>>>>>> "+time.ctime()        
        commands.getoutput("3dcalc -a "+self.options.loc+"/"+self.options.id+".around.vent.blur+orig -b "+frac_mask+" -c "+self.options.loc+"/vent.seed."+self.options.id+"+orig -expr 'step(b)*iszero(step(a))*(1+100*step(c))' -prefix "+self.options.loc+"/"+self.options.id+".around.vent.inv")
        print "<<<<<<<< 3dmerge to cluster around the vent >>>>>>>>>>>>> "+time.ctime()
        commands.getoutput("3dmerge -dxyz=1 -1clust_max 2 1 -prefix "+self.options.loc+"/"+self.options.id+".around.vent.clust "+self.options.loc+"/"+self.options.id+".around.vent.inv+orig")
        print "<<<<<<<< 3dcalc to get masked volume ventricles >>>>>>>>>>>>> "+time.ctime()
        commands.getoutput("3dcalc -datum byte -a "+self.options.loc+"/"+self.options.id+".around.vent.clust+orig -expr '100*step(a-1)' -prefix "+self.options.loc+"/"+self.options.id+".vent.init")
        print "<<<<<<<< 3dmerge to get masked volume ventricles blur >>>>>>>>>>>>> "+time.ctime()
        commands.getoutput("3dmerge -1filter_nzmean 5 -prefix "+self.options.loc+"/"+self.options.id+".vent.init.blur "+self.options.loc+"/"+self.options.id+".vent.init+orig")
        print "<<<<<<<< 3dcalc to grab masked volume ventricles using the blurred >>>>>>>>>>>>> "+time.ctime()   
        commands.getoutput("3dcalc -a "+self.options.loc+"/"+self.options.id+".around.vent.preblur+orig -b "+self.options.loc+"/"+self.options.id+".vent.init.blur+orig -c "+frac_mask+" -expr 'iszero(a)*iszero(b)*step(c)' -prefix "+self.options.loc+"/"+self.options.id+".around.vent.init")
        print "<<<<<<<< 3dmerge to get init volume around ventricles >>>>>>>>>>>>> "+time.ctime()
        commands.getoutput("3dmerge -dxyz=1 -1clust_order 2 1 -prefix "+self.options.loc+"/"+self.options.id+".around.vent.init.clust "+self.options.loc+"/"+self.options.id+".around.vent.init+orig")
        print "<<<<<<<< 3dcalc to get init preclustered volume around ventricles >>>>>>>>>>>>> "+time.ctime()     
        commands.getoutput("3dcalc -a "+self.options.loc+"/"+self.options.id+".around.vent.init.clust+orig -b "+self.options.loc+"/"+self.options.id+".vent.init+orig -c "+self.options.loc+"/"+self.options.id+".around.vent.preblur+orig -d "+frac_mask+" -expr 'iszero(c)*iszero(equals(a,1))*(1+100*step(b))*step(d)' -prefix "+self.options.loc+"/"+self.options.id+".vent.init.preclust")
        print "<<<<<<<< 3dmerge to get initial clustered volume around ventricles >>>>>>>>>>>>> "+time.ctime()
        commands.getoutput("3dmerge -dxyz=1 -1clust_max 2 1 -prefix "+self.options.loc+"/"+self.options.id+".vent.init.clust "+self.options.loc+"/"+self.options.id+".vent.init.preclust+orig")        
        ## finally getting the ventricles mask
        print "<<<<<<<< 3dcalc and 3dfractionize to get ventricles mask >>>>>>>>>>>>> "+time.ctime()
        os.system("3dcalc -datum byte -a "+self.options.loc+"/"+self.options.id+".vent.init.clust+orig -expr 'step(a-1)' -prefix "+self.options.loc+"/volume.VENT."+self.options.id)
        print "<<<<<<<< 3dcalc to get white matter mask >>>>>>>>>>>>> "+time.ctime()
        commands.getoutput("3dcalc -datum byte -a "+self.options.loc+"/volume."+self.options.id+"_seg.nii.gz -expr 'equals(a,3)' -prefix "+self.options.loc+"/volume.WM."+self.options.id)

    def Clipper(self):
        print """Now get quartile clips"""
        auto_mask = self.options.mask
        vent_med1 = commands.getoutput("3dmaskave -median -mask "+self.options.loc+"/volume.VENT."+self.options.id+"+orig "+self.options.vol+" | awk '{print $1}'").split('\n')[2]
        print "vent_med1: "+vent_med1
        vent_med2 = commands.getoutput("3dmaskave -median -mask "+self.options.loc+"/volume.VENT."+self.options.id+"+orig -drange 0 "+vent_med1+" "+self.options.vol+" | awk '{print $1}'").split('\n')[2]
        print "vent_med2: "+vent_med2
        print "now running: 3dcalc -a "+self.options.vol+" -b "+self.options.loc+"/volume.VENT."+self.options.id+"+orig -expr 'step("+vent_med2+"+1-a)*step(b)' -prefix "+self.options.loc+"/mask.VENTclip."+self.options.id+"+orig "+time.ctime()
        commands.getoutput("3dcalc -a "+self.options.vol+" -b "+self.options.loc+"/volume.VENT."+self.options.id+"+orig. -expr 'step("+vent_med2+"+1-a)*step(b)' -prefix "+self.options.loc+"/mask.VENTclip."+self.options.id+"+orig") ##vent_med2+1 is used to prevent restricted ranges from not producing masks
        wm_med1 = commands.getoutput("3dmaskave -median -mask "+self.options.loc+"/volume.WM."+self.options.id+"+orig "+self.options.vol+" | awk '{print $1}'").split('\n')[2] 
        print "wm_med1: "+wm_med1
        wm_med2 = commands.getoutput("3dmaskave -median -mask "+self.options.loc+"/volume.WM."+self.options.id+"+orig -drange "+wm_med1+" 100000 "+self.options.vol+" | awk '{print $1}'").split('\n')[2] 
        print "wm_med2: "+wm_med2 
        commands.getoutput("3dcalc -a "+self.options.vol+" -b "+self.options.loc+"/volume.WM."+self.options.id+"+orig. -expr 'step(a-"+wm_med2+")*step(b)' -prefix "+self.options.loc+"/mask.WMclip."+self.options.id+"+orig") 
        for i in self.options.nn.split('\t')[:]:
            idname = self.options.id+"."+i
            commands.getoutput("3dresample -master "+self.options.path+"/reg.despike.shift."+idname+"+orig. -prefix "+self.options.loc+"/automask."+idname+" -inset "+self.options.mask)
            commands.getoutput("3dfractionize -template "+self.options.loc+"/automask."+idname+"+orig. -input "+self.options.loc+"/mask.VENTclip."+self.options.id+"+orig -prefix "+self.options.loc+"/mask.VENTclip.frac."+idname+"+orig -vote -clip 0.8")#changed from -clip 0.5 (looks better) 
            commands.getoutput("3dfractionize -template "+self.options.loc+"/automask."+idname+"+orig. -input "+self.options.loc+"/mask.WMclip."+self.options.id+"+orig. -prefix "+self.options.loc+"/mask.WMclip.frac."+idname+"+orig -vote -clip 0.9")
            median_vals = idname+"\nwm_med1: "+wm_med1+"\nwm_med2: "+wm_med2+"\nvent_med1: "+vent_med1+"\nvent_med2: "+vent_med2+"\n" #recording wm and vent values 
            file = open(self.options.loc+"/median_vals.1D",'a')
            file.write(median_vals)
            file.close()

    def CheckAndClean(self):
        os.chdir(self.options.loc)
        for i in self.options.nn.split('\t')[:]:
            idname = self.options.id+"."+i
            if exists (self.options.loc+"/mask.VENTclip.frac."+idname+"+orig.BRIK")== False:
                print "\n\n\nERROR: This script has NOT successfully completed. Check subject and file names and data location, etc. Look at intermediary files and determine error source. Error occurred in producing clipped ventricle mask for "+idname+".\n\n\n"
                quit()
            if exists (self.options.loc+"/mask.WMclip.frac."+idname+"+orig.BRIK")== False:
                print "\n\n\nERROR: This script has NOT successfully completed. Look at intermediary files and determine error source. Error occurred in producing clipped white matter mask for "+idname+".\n\n\n"
                quit()
        #if exists (self.options.loc+"/Maskmaker_files")== False:
            #commands.getoutput("mkdir Maskmaker_files")
        #commands.getoutput("mv *vent* tal2volume* *volume.padded."+self.options.id+"* template.resampled* *nii* volume.masked.* "+idname+".automask_frac./Maskmaker_files")
        print "\n\n\nCongratulations! This script (5.Maskmaker) seems to have run successfully! Check the white matter and ventricle masks which have been created by viewing them as overlays "+self.options.vol+", from which they were created. Please check this!\n\n\nIt is not uncommon for this script to be problematic... Check location of vent seed and blurring over ventricles. Good luck!\n\n\n"



mm = MaskMaker()
mm.get_opts()
mm.TLRC_algn()
mm.CalcMaskedVol()
mm.CalcMaskedVolVent()
mm.Clipper()
mm.CheckAndClean()
