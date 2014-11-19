#!/usr/bin/python

import os
from os.path import exists
import commands
import sys
import utils
import time
from optparse import OptionParser


class AlignShiftSpike:
   
	#needs to make directory and needs to do concat on down
	
    def AlignEPICenter(self, epipath, epiname):
	command="""3drefit -xorigin cen -yorigin cen -zorigin cen %s/%s""" %(epipath, epiname)
	os.system(command)


    def Tshift(self, sliceordering, epipath, epiname):
	command="""3dTshift -tpattern %s -ignore 5 -prefix %s/3.Register_Despike/3.ts.shift.%s""" %(sliceordering, epipath, epiname)
	os.system(command)
	

    def Despike(self, epipath, epiname):
	command="""3dDespike -nomask -ignore 5 -prefix %s/3.Register_Despike/3.despike.shift.%s""" %(epipath, epiname)
	os.sytem(command)	

    def Resample(self, epipath, epiname):
	command= """3dresample -master 3.Register_Despike/%s -prefix 3.Register_Despike/3.resampled.%s -inset 3.Register_Despike/3.despike.shift.%s""" %(epipath, epiname, epiname)
	os.system(command)

#I have to see this for myself, I think we can just do a cat or an "&"
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
