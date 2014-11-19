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

class RegSpikePlot:
   
    def get_opts(self):
        desc="""This script registers the functional data for all runs to a pre-selected subbrick from one run (specified in directory pathfile), then despikes the time series and produces a motion plot for each individual condition run and all runs together."""   
        self.usage = "usage: %prog [options]"
        self.parser = OptionParser(description=desc, version="%prog 28.April.2010")
        self.parser.add_option("--subject", dest="ss",
            help="subject identifier-- must match label used for subject directory and HEAD/BRIK files")
        self.parser.add_option("--epi", dest="epi",
            help="functional run to which all runs will be registered (from individual pathfiles)")
        self.parser.add_option("--epibase", dest="epibase",
            help="from individual pathfiles- this is the selected TR from the epi")
        self.parser.add_option("--newnames", dest="nn",
            help="suffix for files for individual runs")
        self.parser.add_option("--path", dest="path",
            help="path to preprocessing data")
        (self.options, args) = self.parser.parse_args()

    def AlignEPICenter(self):
        os.chdir(self.options.path)
	commands.getoutput("mkdir 3.Register_Despike")
	print_out(self, "--- Aligning center of run to center of coordinate space")
        for i in self.options.nn.split('\t')[:]:
            idname = self.options.ss+"."+i
	    print "\033[30m \n \tScan: " + idname
            commands.getoutput("3drefit -xorigin cen -yorigin cen -zorigin cen 1.ConvertMR/1.ts." + idname + "+orig.BRIK")

	#ASK SUSAN:  SHOULDN'T THE LINE ABOVE BE INDENTED WITH THE FOR LOOP?

    def Tshift(self):
        os.chdir(self.options.path)
	print_out(self, "--- Performing slice time correction")
        for i in self.options.nn.split('\t')[:]:
            idname = self.options.ss+"."+i
	    print "\033[30m \n \tScan: " + idname
	    #in next line, "-ignore 5" discards first 5 crummy scans...check to see if increasing DS from 2-7 fixes this.
            commands.getoutput("3dTshift -tpattern altplus -ignore 5 -prefix 3.Register_Despike/3.ts.shift."+idname+" 1.ConvertMR/1.ts."+idname+"+orig.")

    def Despike(self):
        os.chdir(self.options.path) # path to ts
	print_out(self, "--- Despiking Time Series to center of cordinate space to center of coordinate space")	
        for i in self.options.nn.split('\t')[:]:
            idname = self.options.ss+"."+i
	    print "\033[30m \n \tScan: " + idname
            #in next line, "-ignore 5" discards first 5 crummy scans...check to see if increasing DS from 2-7 fixes this.
            commands.getoutput("3dDespike -nomask -ignore 5 -prefix 3.Register_Despike/3.despike.shift."+idname+" 3.Register_Despike/3.ts.shift."+idname+"+orig")

    def Volreg(self):
        os.chdir(self.options.path)
	print_out(self, "--- Resampling the base time series to the space of the EPI being registered to it \n--- Registering EPI runs to base ")	
        for i in self.options.nn.split('\t')[:]:
            idname = self.options.ss+"."+i
####2014-10-06: ESD/ADS: The folowing lines have been commented/amended as all fMRI runs use the same voxel size 3X3X3 and resampling is therefore unnecessary. Also, it was fucking shit up (i.e. clipping off part of the brain (posterior).
	    print "\033[30m \n \tScan: " + idname
	    commands.getoutput("3dresample -master 3.Register_Despike/3."+ self.options.epi+"+orig. -prefix 3.Register_Despike/3.resampled."+idname+" -inset 3.Register_Despike/3.despike.shift."+idname+"+orig.BRIK")

            ##commands.getoutput("3dvolreg -overwrite -twopass -twodup -Fourier -zpad 20 -dfile motion_"+idname+" -base 'resampled."+self.options.epi+"+orig.BRIK["+self.options.epibase+"]' -prefix reg.despike.shift."+idname+" despike.shift."+idname+"+orig.BRIK")
#2014-10-06: ADS: i think the problem is that here we don't use "resampled."+self.options.epi+" as the argument for the base flag print commands.getoutput("3dvolreg -twopass -twodup -Fourier -zpad 20 -dfile motion_"+idname+" -base '"+self.options.epi+"+orig.BRIK["+self.options.epibase+"]' -prefix reg.despike.shift."+idname+" despike.shift."+idname+"+orig.BRIK")
####so i commented out the original 3dvolreg command above and changed the base flag argument in the line below....
	
	    print commands.getoutput("3dvolreg -twopass -twodup -Fourier -zpad 10 -dfile 3.Register_Despike/3.motion_" + idname +  " -base '3.Register_Despike/3."+self.options.epi+"+orig.BRIK["+self.options.epibase+"]' -prefix 3.Register_Despike/3.reg.despike.shift."+idname+" 3.Register_Despike/3.resampled."+idname+"+orig.BRIK")
#            commands.getoutput("rm -rf resampled."+self.options.epi+"*")

    def Concat(self):
        os.chdir(self.options.path)
	print_out(self, "--- Concatenating motion parameters")
        for i in self.options.nn.split('\t')[:]:
            idname = self.options.ss+"."+i
            newmotion=commands.getoutput("1dcat 3.Register_Despike/3.motion_"+idname+"[1-6]")
	    print "\nmotion_"+idname+".1D"
            file = open("3.Register_Despike/3.motion_"+idname+".1D",'w')
            file.write(newmotion)
            file.close
            m0=commands.getoutput("head -n 1 3.Register_Despike/3.motion_"+idname+" | awk '{print $2}'")
            print "m0 = "+`m0`
            newm0=commands.getoutput("1deval -a 3.Register_Despike/3.motion_"+idname+"[1] -expr 'a-"+m0+"'")
            file=open("3.Register_Despike/3.mot_m0."+idname+".1D",'w')
            file.write(newm0)
            file.close
            m1=commands.getoutput("head -n 1 3.Register_Despike/3.motion_"+idname+" | awk '{print $3}'")
            print "m1 = "+`m1`
            newm1=commands.getoutput("1deval -a 3.Register_Despike/3.motion_"+idname+"[2] -expr 'a-"+m1+"'")
            file=open("3.Register_Despike/3.mot_m1."+idname+".1D",'w')
            file.write(newm1)
            file.close
            m2=commands.getoutput("head -n 1 3.Register_Despike/3.motion_"+idname+" | awk '{print $4}'")
            print "m2 = "+`m2`
            newm2=commands.getoutput("1deval -a 3.Register_Despike/3.motion_"+idname+"[3] -expr 'a-"+m2+"'")
            file=open("3.Register_Despike/3.mot_m2."+idname+".1D",'w')
            file.write(newm2)
            file.close
            m3=commands.getoutput("head -n 1 3.Register_Despike/3.motion_"+idname+" | awk '{print $5}'")
            print "m3 = "+`m3`
            newm3=commands.getoutput("1deval -a 3.Register_Despike/3.motion_"+idname+"[4] -expr 'a-"+m3+"'")
            file=open("3.Register_Despike/3.mot_m3."+idname+".1D",'w')
            file.write(newm3)
            file.close
            m4=commands.getoutput("head -n 1 3.Register_Despike/3.motion_"+idname+" | awk '{print $6}'")
            print "m4 = "+`m4`
            newm4=commands.getoutput("1deval -a 3.Register_Despike/3.motion_"+idname+"[5] -expr 'a-"+m4+"'")
            file=open("3.Register_Despike/3.mot_m4."+idname+".1D",'w')
            file.write(newm4)
            file.close
            m5=commands.getoutput("head -n 1 3.Register_Despike/3.motion_"+idname+" | awk '{print $7}'")
            print "m5 = "+`m5`
            newm5=commands.getoutput("1deval -a 3.Register_Despike/3.motion_"+idname+"[6] -expr 'a-"+m5+"'")
            file=open("3.Register_Despike/3.mot_m5."+idname+".1D",'w')
            file.write(newm5)
            file.close

    def PlotAll(self):
        os.chdir(self.options.path)
        length=range(len(self.options.nn.split('\t')))
        for i in self.options.nn.split('\t')[:]:
            idname = self.options.ss+"."+i
            file=open("tmp_plot",'a')
            file.write("3.Register_Despike/3.motion_"+idname+" ")
            file.close
	print_out(self, "--- Creating Motion Plot (all conditions)")        
        allmotion=open("tmp_plot").read()        
        print commands.getoutput("cat "+allmotion+" > 3.Register_Despike/3.motion_"+self.options.ss+"_all")
        commands.getoutput("1dplot -ps -volreg -title '''Motion for All Runs of "+self.options.ss+"''' -xlabel 'Time Point' '3.Register_Despike/3.motion_"+self.options.ss+"_all[1-6]' > 3.Register_Despike/3.motion_"+self.options.ss+".ps")
        commands.getoutput("rm -rf tmp_plot")
   
    def MotionPlot(self):
        os.chdir(self.options.path)
	print_out(self, "--- Creating Motion Plots (Individual Conditions)") 
        for i in self.options.nn.split('\t')[:]:
            idname = self.options.ss+"."+i
	    print "\033[30m \n \t" + idname
            commands.getoutput("1d_tool.py -set_nruns 1 -overwrite -infile 3.Register_Despike/3.motion_"+idname+".1D -derivative -write 3.Register_Despike/3.derivative."+idname+".1D")
            tr = commands.getoutput("3dinfo 3.Register_Despike/3.reg.despike.shift."+idname+"+orig | grep '''Time step = ''' | awk -F'''Time step = ''' '{print $2}' | awk -F's' '{print $1}'").split('\n')[-1]
            print "\tTR = "+`tr`

            commands.getoutput("1d_tool.py -set_nruns 1 -overwrite -infile 3.Register_Despike/3.motion_"+idname+".1D -set_tr "+`tr`+" -censor_prev_TR -censor_motion 1 "+ "3.Register_Despike/3."+idname)
	    
            commands.getoutput("1dcat '3.Register_Despike/3.mot_m0."+idname+".1D' '3.Register_Despike/3.mot_m1."+idname+".1D' '3.Register_Despike/3.mot_m2."+idname+".1D' '3.Register_Despike/3.mot_m3."+idname+".1D' '3.Register_Despike/3.mot_m4."+idname+".1D' '3.Register_Despike/3.mot_m5."+idname+".1D' > 3.Register_Despike/3.motion_"+idname+".zerostart.1D")
            commands.getoutput("rm -f 3.Register_Despike/3.mot_m0 3.Register_Despike/3.mot_m1 3.Register_Despike/3.mot_m2 3.Register_Despike/3.mot_m3 3.Register_Despike/3.mot_m4 3.Register_Despike/3.mot_m5")
            commands.getoutput("1dplot -ps -volreg -title '''Motion for "+i+" of "+self.options.ss+"''' -ylabel '(Initial shift has been subtracted)' -xlabel 'Time Point' 3.Register_Despike/3.motion_"+idname+".zerostart.1D 3.Register_Despike/3."+idname+"_censor.1D > 3.Register_Despike/3.motion."+idname+".ps")
	    commands.getoutput("rm -rf 3.Register_Despike/3.mot_*"+idname+".1D")

    def CheckAndClean(self):
        os.chdir(self.options.path)
        for i in self.options.nn.split('\t')[:]:
            idname = self.options.ss+"."+i
            if exists (self.options.path+"/3.Register_Despike/3.ts.shift."+idname+"+orig.BRIK") == False:
                print "\n" + Tcolor.ERROR + "ERROR: This script has NOT successfully completed. Check subject and file names and data location, etc. Slice timing correction (and all following steps) failed. Error occurred on "+idname+"." + Tcolor.END + "\n"
                quit()
            if exists (self.options.path+"/3.Register_Despike/3.despike.shift."+idname+"+orig.BRIK") == False:
                print "\n" + Tcolor.ERROR + "ERROR: This script has NOT successfully completed. Check subject and file names and data location, etc. Despiking (and all following steps) failed. Error occurred on "+idname+"." + Tcolor.END + "\n"
                quit()
            if exists (self.options.path+"/3.Register_Despike/3.reg.despike.shift."+idname+"+orig.BRIK")== False:
                print "\n" + Tcolor.ERROR + "ERROR: This script has NOT successfully completed. Check subject and file names and data location, etc. Motion registration (and all following steps) failed. Error occurred on "+idname+"." + Tcolor.END + "\n"
                quit()
            if exists (self.options.path+"/3.Register_Despike/3.motion."+idname+".ps")== False:
                print "\n" + Tcolor.ERROR + "ERROR: This script has NOT successfully completed. Check subject and file names and data location, etc. Motion plots have not been created. Error occurred on "+idname+". Hopefully this doesn't mean there's a problem with your motion regressors...." + Tcolor.END + "\n"
                quit()  
        print "\n" + Tcolor.SUCCESS + "Congratulations! This script (3.registerdespikeplot) was a success! The functional MR runs have been registered to the selected reference scan (resampled to allow for differences in voxel size) and reference TR ('"+self.options.epi+"+orig.BRIK["+self.options.epibase+"]). Motion plots and censor files have also been created and can (should!) be viewed." + Tcolor.END + "\n"

                                  
rsp = RegSpikePlot()
rsp.get_opts()
rsp.AlignEPICenter()
rsp.Tshift()
rsp.Despike()
rsp.Volreg()
rsp.Concat()
rsp.PlotAll()
rsp.MotionPlot()
rsp.CheckAndClean()
