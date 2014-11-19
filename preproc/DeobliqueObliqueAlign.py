class DeobliqueObliqueAlign:

#why are we running this?
#Instead of running two commanads
#commands.getoutput("mkdir 4.Deoblique_Oblique_Align")
    def Deoblique(self, anatname, path):
	
		os.chdir(path) # path to volume
        command="""3dWarp -deoblique -zpad 5 -prefix 4.Deoblique_Oblique_Align/4.volume.deobl.stripped.%s -gridset 2.Skullstrip_Zeropad/2.volume.padded.stripped.%s  2.Skullstrip_Zeropad/2.volume.padded.stripped.%s""" %(anatname, anatname, anatname) # presumes volume has been padded
		os.system(command)
        command="""3dWarp -deoblique -zpad 5 -prefix 4.Deoblique_Oblique_Align/4.volume.deobl.skull.%s -gridset 2.Skullstrip_Zeropad/2.volume.padded.skull.%s 2.Skullstrip_Zeropad/2.volume.padded.skull.%s"""%(anatname, anatname, anatname) # presumes volume has been padded
		os.system(command)

    def Oblique(self,path, epiname,anatname):
        print_out(self, "--- Oblique started") #what is this doing? to what?
        os.chdir(path) # path to volume	
        command="""3dWarp -oblique_parent 3.Register_Despike/3.%s -zpad 5 -prefix 4.Deoblique_Oblique_Align/4.volume.oblepi.stripped.%s -gridset 4.Deoblique_Oblique_Align/4.volume.deobl.stripped.%s 4.Deoblique_Oblique_Align/4.volume.deobl.stripped.%s""" %(epiname,anatname, anatname, anatname)
		os.system(command)
        command="""3dWarp -oblique_parent 3.Register_Despike/3.%s -zpad 5 -prefix 4.Deoblique_Oblique_Align/4.volume.oblepi.skull.%s -gridset 4.Deoblique_Oblique_Align/4.volume.deobl.skull.%s 4.Deoblique_Oblique_Align/4.volume.deobl.skull.%s""" %(epiname,anatname, anatname, anatname)
		os.system(command)

    def Align(self,path, epiname,epibasename, anatname):
		print_out(self, "--- Alignment started")
        os.chdir(path) # path to volume
		command="""align_epi_anat.py -ex_mode quiet -anat2epi -deoblique off -volreg off -anat 4.Deoblique_Oblique_Align/4.volume.oblepi.stripped.%s -cost lpc -anat_has_skull no -AddEdge -overwrite -epi 3.Register_Despike/3.%s -epi_base %s -epi_strip 3dAutomask -master_anat SOURCE -cmass cmass -big_move -child_anat 4.Deoblique_Oblique_Align/4.volume.oblepi.skull.%s""" %(anatname, epiname, epibasename, anatname)
		os.system(command)
		#commands.getoutput("mv 4.* 4.Deoblique_Oblique_Align") 

    def Automask(self, path, anatname):
		print_out(self, "--- Making automask")
        os.chdir(path)# path to volume
		command="""3dAutomask -prefix 4.Deoblique_Oblique_Align/4.anat_automask.%s -dilate 0 4.Deoblique_Oblique_Align/4.volume.oblepi.stripped.%s_al+orig.""" %(anatname, anatname)
		os.system(command)

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
