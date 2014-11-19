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

class DeobliqueObliqueAlign:

    def get_opts(self):
        desc="""This script deobliques and then obliques the volume. Check both to make sure no clipping occurs (if so more zeropadding is needed)."""
        self.usage = "usage: %prog [options]"
        self.parser = OptionParser(description=desc, version="%prog 28.April.2010")
        self.parser.add_option("--subject", dest="ss",
            help="subject identifier-- must match label used for subject directory and HEAD/BRIK files")
        self.parser.add_option("--epi", dest="epi",
            help="from individual pathfiles")
        self.parser.add_option("--epibase", dest="epibase",
            help="from individual pathfiles")
        self.parser.add_option("--newnames", dest="nn",
            help="suffix for files for individual runs")
        self.parser.add_option("--path", dest="path",
            help="path to preprocessing directory")
        (self.options, args) = self.parser.parse_args()

    def Deoblique(self):
	print_out(self, "--- Deoblique started")
        os.chdir(self.options.path) # path to volume
	commands.getoutput("mkdir 4.Deoblique_Oblique_Align")
        commands.getoutput("3dWarp -deoblique -zpad 5 -prefix 4.Deoblique_Oblique_Align/4.volume.deobl.stripped."+self.options.ss+" -gridset 2.Skullstrip_Zeropad/2.volume.padded.stripped."+self.options.ss+"+orig 2.Skullstrip_Zeropad/2.volume.padded.stripped."+self.options.ss+"+orig") # presumes volume has been padded
        commands.getoutput("3dWarp -deoblique -zpad 5 -prefix 4.Deoblique_Oblique_Align/4.volume.deobl.skull."+self.options.ss+" -gridset 2.Skullstrip_Zeropad/2.volume.padded.skull."+self.options.ss+"+orig 2.Skullstrip_Zeropad/2.volume.padded.skull."+self.options.ss+"+orig") # presumes volume has been padded

    def Oblique(self):
        print_out(self, "--- Oblique started")
        os.chdir(self.options.path) # path to volume	
        commands.getoutput("3dWarp -oblique_parent " + "3.Register_Despike/3." + self.options.epi+"+orig.BRIK -zpad 5 -prefix 4.Deoblique_Oblique_Align/4.volume.oblepi.stripped."+self.options.ss+" -gridset 4.Deoblique_Oblique_Align/4.volume.deobl.stripped."+self.options.ss+"+orig 4.Deoblique_Oblique_Align/4.volume.deobl.stripped."+self.options.ss+"+orig")
        commands.getoutput("3dWarp -oblique_parent " + "3.Register_Despike/3." + self.options.epi+"+orig.BRIK -zpad 5 -prefix 4.Deoblique_Oblique_Align/4.volume.oblepi.skull."+self.options.ss+" -gridset 4.Deoblique_Oblique_Align/4.volume.deobl.skull."+self.options.ss+"+orig 4.Deoblique_Oblique_Align/4.volume.deobl.skull."+self.options.ss+"+orig")

    def Align(self):
	print_out(self, "--- Alighnment started")
        os.chdir(self.options.path)# path to volume
	commands.getoutput("align_epi_anat.py -ex_mode quiet -anat2epi -deoblique off -volreg off -anat 4.Deoblique_Oblique_Align/4.volume.oblepi.stripped."+self.options.ss+"+orig.BRIK -cost lpc -anat_has_skull no -AddEdge -overwrite -epi 3.Register_Despike/3."+self.options.epi+"+orig.BRIK -epi_base "+self.options.epibase+" -epi_strip 3dAutomask -master_anat SOURCE -cmass cmass -big_move -child_anat 4.Deoblique_Oblique_Align/4.volume.oblepi.skull."+self.options.ss+"+orig.BRIK")
	commands.getoutput("mv 4.* 4.Deoblique_Oblique_Align") 

    def Automask(self):
	print_out(self, "--- Making automask")
        os.chdir(self.options.path)# path to volume
	commands.getoutput("3dAutomask -prefix 4.Deoblique_Oblique_Align/4.anat_automask."+self.options.ss+" -dilate 0 4.Deoblique_Oblique_Align/4.volume.oblepi.stripped."+self.options.ss+"_al+orig.")
#        for di in range(1,5):
#        commands.getoutput("3dAutomask -prefix anat_automask."+self.options.ss+" -dilate 4 volume.oblepi.stripped."+self.options.ss+"_al+orig.")
#        for i in self.options.nn.split('\t')[:]:
#            idname = self.options.ss+"."+i
#            commands.getoutput("3dresample -master reg.despike.shift."+idname+"+orig. -prefix automask."+idname+" -inset anat_automask."+self.options.ss+"+orig.")           

    def CheckAndClean(self):
        os.chdir(self.options.path)
        if exists (self.options.path+"/4.Deoblique_Oblique_Align/4.volume.deobl.stripped."+self.options.ss+"+orig.BRIK")== False:
            print "\n" + Tcolor.ERROR + "ERROR: This script has NOT successfully completed. Check subject and file names and data location, etc. Error occurred in deobliquing "+self.options.ss+"." + Tcolor.END + "\n"
            quit()
        if exists (self.options.path+"/4.Deoblique_Oblique_Align/4.volume.oblepi.stripped."+self.options.ss+"+orig.BRIK")== False:
            print "\n" + Tcolor.ERROR + "ERROR: This script has NOT successfully completed. Check subject and file names and data location, etc. Error occurred in re-obliquing "+self.options.ss+"." + Tcolor.END + "\n"
            quit()
        if exists (self.options.path+"/4.Deoblique_Oblique_Align/4.volume.oblepi.stripped."+self.options.ss+"_al+orig.BRIK")== False:
            print "\n" + Tcolor.ERROR + "ERROR: This script has NOT successfully completed. Check subject and file names and data location, etc. Error occurred in aligning "+self.options.ss+" to "+self.options.epi+"." + Tcolor.END + "\n"
            quit()
#        for di in range(1,5):
        if exists (self.options.path+"/4.Deoblique_Oblique_Align/4.anat_automask."+self.options.ss+"+orig.BRIK")== False:
            print "\n" + Tcolor.ERROR + "ERROR: This script has NOT successfully completed. Check subject and file names and data location, etc. Error occurred in constructing automask for "+self.options.ss+"." + Tcolor.END + "\n"
            quit()
#        for i in self.options.nn.split('\t')[:]:
#            idname= self.options.ss+"."+i
#            if exists (self.options.path+"/automask."+idname+"+orig.BRIK")==False:
#                print "\n\nERROR: This script has NOT successfully completed. Check subject and file names and data location, etc. Error occurred in constructing automask for "+idname+".\n\n"
#                quit()
	commands.getoutput("mv AddEdge 4.Deoblique_Oblique_Align/4.AddEdge")        
	#os.chdir("AddEdge")
        #commands.getoutput("3drename volume.oblepi.stripped."+self.options.ss+"_al_rs_ec volumewithedge."+self.options.ss+"_al") --may want to add these back in to simplify addedge comparison-- but check which are best for visualizing
        #commands.getoutput("3drename "+self.options.epi+"_ns_ec epiwithedge."+self.options.epi)
        #commands.getoutput("mv "+self.options.path+"/AddEdge/ts.* "+self.options.path)
        #commands.getoutput("mv "+self.options.path+"/AddEdge/volume.oblepi.* "+self.options.path)
	print "\n" + Tcolor.SUCCESS + "Congratulations! This script (4.deoblique_oblique_align) was a success! The anatomical volume has been deobliqued from the angle at which it was acquired (to plumb) and then re-obliqued and aligned to the functional run selected to serve as reference ("+self.options.epi+").\n\nYou can (and should!) check the accuracy of the alignment. This can be done by overlaying the _al anatomical on the fmri and vice versa. You should also check the anatomical and functional scans in "+self.options.path+"/AddEdge. If it doesn't look good.... well who know? Perhaps another cost function (this script uses Local Pearson Correlation) may be better for this data set? Tweak away! Good luck!" + Tcolor.END + "\n"


doa=DeobliqueObliqueAlign()
doa.get_opts()
doa.Deoblique()
doa.Oblique()
doa.Align()
doa.Automask()
doa.CheckAndClean()



#    def Automask(self):
#        print "\nMaking automask for "+self.options.ss+"---"+time.ctime()
#        os.chdir(self.options.path)# path to volume
#        for di in range(1,5):
#            commands.getoutput("3dAutomask -prefix anat_automask."+self.options.ss+"."+`di`+" -dilate "+`di`+" volume.oblepi.stripped."+self.options.ss+"_al+orig.")
#            for i in self.options.nn.split('\t')[:]:
#                idname = self.options.ss+"."+i
#                commands.getoutput("3dresample -master reg.despike.shift."+idname+"+orig. -prefix automask."+idname+"."+`di`+" -inset anat_automask."+self.options.ss+"."+`di`+"+orig.") 

#        for di in range(1,5):
#            if exists (self.options.path+"/anat_automask."+self.options.ss+"."+`di`+"+orig.BRIK")== False:
#                print "\n\nERROR: This script has NOT successfully completed. Check subject and file names and data location, etc. Error occurred in constructing automask for "+self.options.ss+".\n\n"
#                quit()
#            for i in self.options.nn.split('\t')[:]:
#                idname= self.options.ss+"."+i
#                if exists (self.options.path+"/automask."+idname+"."+`di`+"+orig.BRIK")==False:
#                    print "\n\nERROR: This script has NOT successfully completed. Check subject and file names and data location, etc. Error occurred in constructing automask for "+idname+".\n\n"
#                    quit()
