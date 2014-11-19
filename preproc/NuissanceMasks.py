#!/usr/bin/python
#NOTE: CHECK BLUR in line 78 and in "+ss+".around.vent.blur!!!!

import os
from os.path import exists
import commands
import sys
import utils
import time
from optparse import OptionParser


class NuissanceMasks:

             
    def TLRC_algn(self, subid, subid, path, anatname):
      
		print_out(self, "--- Zeropad and get a tal align coordinate vector") #why are we zeropadding?
        command="""3dZeropad -I 10 -S 10 -A 10 -P 10 -L 10 -R 10 -prefix %s/volume.padded.%s %s""" %(path, subid, anatname) ## create extra padding in the volume space### 
		os.system(command)
		commands="""3dresample -inset /soft/afni/TT_N27+tlrc -master %s/volume.padded.%s+orig -prefix %s/template.resampled.%s+orig""" %(path, subid, path,  subid) ## resample to tal
        os.system(command)
		
        print_out(self, "<<<<<<<<<<<<<<<<<<<<<< 3dWarpDrive >>>>>>>>>>>>>>>>>>>>>>>>>>>" )
        commands="""3dWarpDrive -affine_general -cubic -input %s/volume.padded.%s+orig -prefix %s/tal.volume.padded.%s -base %s/template.resampled.%s+orig -1Dmatrix_save %s/tal2volume.%s.1D""" %(path, anatname, path, anatname, path, anatname, path, anatname) ## using this to get the tal transform matrix 'tal2volume.*.1D'
		os.system(command)
		
    def CalcMaskedVol(self, path, subid, anatname, automask, fracmask ):
	
        print_out(self,"Get masked data from the volume for ventricles and white matter.")
        #auto_mask = self.options.mask
        #frac_mask = self.options.loc+"/"+self.options.id+".automask_frac+orig"
        print_out(self, "<<<<<<<<<<<<<<<<<<<<<< 3dfractionize >>>>>>>>>>>>>>>>>>>>>>>>>>> ")
        command="""3dfractionize -template %s -input %s -preserve -clip 0.2 -prefix %s""" %(anatname, automask, fracmask)
        os.system(command)
		
        #there may be need to source FSL but why not include in people's .bashrc... maybe a frontalis issue? 
		#commands.getoutput("source /soft/fsl/etc/fslconf/fsl.sh") ## source FSL
        
		print_out(self, "<<<<<<<< 3dcalc to get masked volume from the fractionized and boxed volume >>>>>>>>>>>>>" )
        command="""3dcalc -a %s -b %s -expr 'step(a)*b' -prefix %path/volume.masked.%s+orig""" %(fracmask, anatname, path, subid) ## mask sample volume
        os.system(command)
		
		print_out(self,"<<<<<<<< 3dAFNItoNIFTI to get nii volume >>>>>>>>>>>>> ") #this should be in the utils.py
        command="""3dAFNItoNIFTI -prefix %s/volume.masked.%s.nii %s/volume.masked.%s+orig""" %(path, subid, path, subid) ## convert to nii
		os.system(command)
		
        print_out(self, "<<<<<<<< running 'fast' to get segmentation (not actually fast-- kinda slow) >>>>>>>>>>>>> ")
        os.system(". ${FSLDIR}/etc/fslconf/fsl.sh")
		command="""fast -o %s/volume.%s %s/volume.masked.%s.nii""" %(path, subid, path,subid)  ## segments volume
        os.system(command)
		
		print_out(self,"<<<<<<<< 3dcalc to get mask around ventricles from nii segmentation >>>>>>>>>>>>> ")
        command="""3dcalc -a %s/volume.%s_seg.nii.gz -b %s -expr '100*step(b)*iszero(amongst(a,1,0))' -prefix %s/%s.around.vent.preblur"""  %(path, subid, fracmask, path, subid) ## mask around the ventricle
		os.system(command)
        
		print_out(self, "<<<<<<<< 3dmerge to blur aroundVent >>>>>>>>>>>>> ")
        ## the blur is not 'one size fits all', i.e., may have to change 1filter_nzmean value so it doesn't blur over the ventricle completely.  
		##start with '-1filter_nzmean 3' per Mike Andric but using 2 here.... seems to work
        command="""3dmerge -1filter_nzmean 2 -prefix %s/%s.around.vent.blur %s/%s.around.vent.preblur+orig""" %s(path,subid,path,subid) ## blur around the ventricle
		os.system(command)
		
    def CalcMaskedVolVent(self, subid, path, anatname, fracmaks):
	
        print_out(self, "<<<<<<<< Generating the Ventricles seed 'Vseed' >>>>>>>>>>>>> ")
        Vseed = "-8 13 19\n8 13 19\n" ## creating seed in areas of ventricle - tal coords (may not be best coords for all- something to check if output is not existent or ventricles incluse area outside skull)
        file = open(path+"/vent.seed.tal.1D",'w')
        file.write(Vseed)
        file.close()
		
        #auto_mask = self.options.mask
        #frac_mask = self.options.loc+"/"+self.options.id+".automask_frac+orig"
		
        print_out(self,"<<<<<<<< Warping from talairach using 'Vecwarp'  >>>>>>>>>>>>> ")
		command="""Vecwarp -matvec %s/tal2volume.%s.1D -forward -input %s/vent.seed.tal.1D -output %s/vent.seed.%s.1D""" %(path, subid, path, path, subid)
        os.system(command)
		
        print_out(self,"<<<<<<<< 3dUndump to get ventricles seed >>>>>>>>>>>>> ")
        command="""3dUndump -xyz -orient RAI -prefix %s/%s" -master %s -srad 8 %s/vent.seed.%s.1D""" %(path, subid, fracmask, path, subid)
		os.system(command)
		
        print_out(self,"<<<<<<<< 3dcalc to get inverted vent volume seed  >>>>>>>>>>>>> ")        
        command="""3dcalc -a %s/%s.around.vent.blur+orig -b %s -c %s/vent.seed.%s+orig -expr 'step(b)*iszero(step(a))*(1+100*step(c))' -prefix %s/%s.around.vent.inv""" %(path, subid, fracmask,path, subid, path, subid)
		
        print_out(self, "<<<<<<<< 3dmerge to cluster around the vent >>>>>>>>>>>>> ")
        command="""3dmerge -dxyz=1 -1clust_max 2 1 -prefix %s/%s.around.vent.clust %s/%s.around.vent.inv+orig""" %(path,subid, path, subid)
		os.system(command)
		
        print_out(self, "<<<<<<<< 3dcalc to get masked volume ventricles >>>>>>>>>>>>> ")
        command="""3dcalc -datum byte -a %s/%s.around.vent.clust+orig -expr '100*step(a-1)' -prefix %s/%s.vent.init""" %(path, subid, path, subid)
		os.system(command)
		
        print_out(self, "<<<<<<<< 3dmerge to get masked volume ventricles blur >>>>>>>>>>>>> ")
        command="""3dmerge -1filter_nzmean 5 -prefix %s/%s.vent.init.blur %s/%s.vent.init+orig""" %(path, subid, path, subid)
		os.system(command)
		
        print_out(self, "<<<<<<<< 3dcalc to grab masked volume ventricles using the blurred >>>>>>>>>>>>> ")   
        command="""3dcalc -a %s/%s".around.vent.preblur+orig -b %s/%s.vent.init.blur+orig -c %s -expr 'iszero(a)*iszero(b)*step(c)' -prefix %s/%s.around.vent.init""" %(path, subid, path, subid, fracmask, path, subid)
		os.system(command)
		
        print_out(self, "<<<<<<<< 3dmerge to get init volume around ventricles >>>>>>>>>>>>> ")
        command="""3dmerge -dxyz=1 -1clust_order 2 1 -prefix %s/%s.around.vent.init.clust %s/%s.around.vent.init+orig""" %(path, subid, path, subid)
		os.system(command)
		
        print_out(self, "<<<<<<<< 3dcalc to get init preclustered volume around ventricles >>>>>>>>>>>>> ")     
        command="""3dcalc -a %s/%s.around.vent.init.clust+orig -b %s/%s.vent.init+orig -c %s/%s.around.vent.preblur+orig -d %s -expr 'iszero(c)*iszero(equals(a,1))*(1+100*step(b))*step(d)' -prefix %s/%s.vent.init.preclust""" %(path, anatname,  path, anatname, fracmask, path, anatname)
		os.system(command)
		
		
        print_out(self, "<<<<<<<< 3dmerge to get initial clustered volume around ventricles >>>>>>>>>>>>> ")
        command="""3dmerge -dxyz=1 -1clust_max 2 1 -prefix %s/%s.vent.init.clust %s/%s.vent.init.preclust+orig""" %(path, subid,path, subid)
	    os.system(command)
	
        ## finally getting the ventricles mask
        print_out(self,"<<<<<<<< 3dcalc and 3dfractionize to get ventricles mask >>>>>>>>>>>>> ")
        command= """3dcalc -datum byte -a %s/%s.vent.init.clust+orig -expr 'step(a-1)' -prefix %s/volume.VENT.%s""" %(path, subid,path, subid)
		os.system(command)
		
		
        print_out(self, "<<<<<<<< 3dcalc to get white matter mask >>>>>>>>>>>>> ")
        command="""3dcalc -datum byte -a %s/volume.%s_seg.nii.gz -expr 'equals(a,3)' -prefix %s/volume.WM.%s""" %(path, subid, path, subid)
		os.system(command)
		

		
		
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
