#!/usr/bin/python

#Script for processing resting data according to Hallquist et al., Neuroimage, 2013
#JRC 11-18-13

#Steps:
# (Note: assume motion correction, slice timing correction, alignment already done and nuisance ROIs already made)
# 1. Smooth (optional)
# 2. Create all nuisance regressors and their derivatives--6 motion parameters, WM, CSF, and optionally global signal, and create one big matrix including all of them
# 3. Bandpass filter EPI data and remove all regressors/derivatives, this includes despiking
# 4. Scrub (optional)
# 5. Correlate (optional)
# NB: If want to concatenate across blocks or conditions, do in a separate script since that will be study-specific

# Usage: python funccon_proc.py project atlas subblocks dosmooth dofilter doscrub docorr
### subblocks = number of subblocks to split each block into--user's responsibility to choose a number of subblocks that divides evenly into the total block length
### dosmooth = smooth alreadysmooth nosmooth
### dofilter = filter alreadyfilter nofilter
### doscrub = scrub noscrub both
### docorr = corr nocorr

# This script requires a params file with columns that are: subject number, num conditions, num blocks per condition, list of conditions

######### IMPORT LIBRARIES #########
import os
import sys
import subprocess
import numpy as np

######### FUNCTIONS #########
def remove_previous_file(filename):
    """ Function that forcibly removes previous versions of a file.
    """

    command = 'rm -f %s' %filename
    os.system(command)
    
def smooth(epidir,epi,epi_smooth,fwhm):
    ''' Smooth EPI data '''

    remove_previous_file('%s%s' %(epidir,epi_smooth)) #added, otherwise doesn't resave
    
    # Note that in fox correlation script smoothing is equivalent to 3dmerge -1blur_fwhm; 3dBlurToFWHM better matches across subs so adopting it here
    #command = """3dBlurToFWHM -input %s/%s -prefix %s/%s -FWHM %s"""
    #%(epidir,epi,epidir,epi_smooth,fwhm) ### REMOVING '/'
    command = """3dBlurToFWHM -input %s%s -prefix %s%s -FWHM %s""" %(epidir,epi,epidir,epi_smooth,fwhm)
    os.system(command)

def get_nuisance_data(s,c,b,mot_gsr_epi_dir,epi_name,nuis_dir,regoutname,GSR=False):
    ''' Get mean and derivative of all motion parameters, white matter/CSF nuisance ROIs, and, optionally, whole brain signal '''

    if GSR==False:
        roi_type = ['mot','white','ventricles']
    else:
        roi_type = ['mot','white','ventricles','gsr']
    reg_type=['deriv','demean']

    for count,r in enumerate(roi_type):

        #define file names
        if r=='mot':
            epi_dir=mot_gsr_epi_dir
            roi_name='dfile.r%02d.1D' %(b)
        elif r=='white' or r=='ventricles':      
            roi_dir=nuis_dir
            if project == 'MegaRest.TMS':
                #need special condition here, not ideal
                #because nuisance ROIs based on baseline sess
                roi_name = '%s_%s_Base.nii' %(r,s[:3])
            else:
                roi_name='%s_%s.nii' %(r,s)
            epi_dir=mot_gsr_epi_dir
        elif r=='gsr':
            roi_dir=mot_gsr_epi_dir
            roi_name='%s-Mask-ALLSESSRUNS.nii' %(s)
            epi_dir=mot_gsr_epi_dir

        #note: will overwrite previous versions of these files (CG)
        if r=='mot': # Just take the motion file output
            # Demean
            command = '1d_tool.py -infile %s/%s -demean -overwrite -write %s/%s_demean.r%02d.1D' %(epi_dir,roi_name,epi_dir,r,b)
            os.system(command)

            # Derivative
            command = '1d_tool.py -infile %s/%s -derivative -demean -overwrite  -write %s/%s_deriv.r%02d.1D' %(epi_dir,roi_name,epi_dir,r,b)
            os.system(command)
            
        else: # Take the average intensity within the ROI first
            # Demean
            command = """3dmaskave -quiet -mask %s/%s %s/%s | 1d_tool.py -infile - -demean -overwrite -write %s/%s_demean.r%02d.1D""" %(roi_dir,roi_name,epi_dir,epi_name,epi_dir,r,b)
            os.system(command)
            
            # Derivative
            command = """3dmaskave -quiet -mask %s/%s %s/%s | 1d_tool.py  -infile - -derivative -demean -overwrite -write %s/%s_deriv.r%02d.1D""" %(roi_dir,roi_name,epi_dir,epi_name,epi_dir,r,b)
            os.system(command)
            
            
        # Now load individual regressor files so they can be concatenated
        command = """%s_demean=np.loadtxt('%s/%s_demean.r%02d.1D')""" %(r,epi_dir,r,b)
        exec(command)
        command = """%s_deriv=np.loadtxt('%s/%s_deriv.r%02d.1D')""" %(r,epi_dir,r,b)
        exec(command)
        
        ### Assuming order of regressors does not matter (demeaned and derivatives interspersed)
        if count==0: # First nuisance regressor so start the nuis_all variable
            if r=='mot': # Only regressor with > 1 column
                command = """nuis_all=np.hstack((%s_demean,%s_deriv))""" %(r,r)
            else: # For regresors with 1 column, need to add a dimension
                command = """nuis_all=np.hstack((%s_demean[:,np.newaxis],%s_deriv[:,np.newaxis]))""" %(r,r)
            exec(command)

        else: # Subsequent nuisance variables add to existing nuis_all
            if r=='mot': # Only regressor with > 1 column
                command = """nuis_all=np.hstack((nuis_all,%s_demean,%s_deriv))""" %(r,r)
            else: # For regresors with 1 column, need to add a dimension
                command = """nuis_all=np.hstack((nuis_all,%s_demean[:,np.newaxis],%s_deriv[:,np.newaxis]))""" %(r,r)
            exec(command)

    # Save variable of all regressors
    regoutfile='%s%s' %(epi_dir,regoutname)
    np.savetxt(regoutfile,nuis_all,fmt='%.6f',delimiter='\t',newline='\n')

def do_bandpass_filtering(epidir,epi,epi_filt,regoutname,filt):
    ''' Do bandpass filtering, regression of nuisance variables, and despiking and output filtered file '''

    remove_previous_file('%s%s' %(epidir,epi_filt)) #added, otherwise doesn't resave

    command = """3dBandpass -despike -ort %s%s -band %.4f %.4f -prefix %s%s -input %s%s""" %(epidir,regoutname,filt[0],filt[1],epidir,epi_filt,epidir,epi)
    os.system(command)

def scrub_data():
    ''' Scrub data to only include low-motion volumes '''
    ##### Holding off on this until determine the best way to execute scrubbing #####
    ##### Note that have to incorporate an option to do both scrubbed and unscrubbed correlations in run_correlations function if doscrub=='both' #####

def run_correlations(maskdir,maskdir_masked,atlas,roilist,maskedroi_prop,mindist,s,epidir,epi_filt,timeslen,start,stop,roi2roidir,timesdir,corroutname,timesoutname,project):
    ''' Run roi2roi correlations '''

    # First define and load anatomical distance mask if atlas is dosenbach_all
    # or power
    nnodes=len(roilist)
    if atlas=='dosenbach_all' or atlas == 'power':
        #maskfilename='%sanatomdist_masks/%s_%s_mindist%s.txt' %(roi2roidir,s,atlas,mindist)
        maskfilename = '%s%s_dist_over%s.txt' %(maskdir,atlas,mindist)
        
        if os.path.exists(maskfilename):
            mask=np.loadtxt(maskfilename)
        else:
            print '.... No distance file found... Assuming all connections can stay, regardless of distance'
            mask=np.ones((nnodes,nnodes)) # If doesn't exist, assume all connections can stay
    else:
        mask=np.ones((nnodes,nnodes)) # Also make a mask of 1s for other atlases to make multiplying easier below
    
        
    for count,r in enumerate(roilist):
        if atlas == 'aal':
            rstr='r%03d' %(r)
        elif atlas == 'Greicius_func' or atlas == 'dosenbach_all' or atlas == 'power':
            rstr='%03d' %(r)
        elif atlas == 'dosenbach' or atlas == 'HO':
            rstr='%02d' %(r)

        #get ROI file names. Clunky to have this if/else here. Should fix.
        if project == 'despolarity':
            roiname =  atlas + "_" + rstr + "_" + s + "r.nii"
            roifile = maskdir + roiname
            roifile_masked = maskdir_masked + roiname
        elif project == 'MegaRest.TMS':
            roiname = '%s_%s_%s' %(atlas,rstr,s[:3])
            roifile = maskdir + roiname + 'r.nii'
            roifile_masked = maskdir_masked + roiname + '.nii'
    

        # Check to make sure enough of masked ROI exists before getting timeseries

        # Full ROI        
        command = """fslstats %s -V""" %(roifile)
        p=subprocess.Popen(command,stdout=subprocess.PIPE,shell=True)
        (tmpout,err) = p.communicate()
        out=tmpout.split()
        roivol=float(out[0]) # Make float so can divide and not output an integer
        
        
        # Masked ROI        
        command = """fslstats %s -V""" %(roifile_masked)
        p=subprocess.Popen(command,stdout=subprocess.PIPE,shell=True)
        (tmpout,err) = p.communicate()
        out=tmpout.split()
        roivol_masked=float(out[0]) # Make float so can divide and not output an integer

        # Check proportion of voxels in masked ROI
        if roivol_masked/roivol < maskedroi_prop:
            # Not enough of ROI left after masking; create a timeseries of nans
            print 'ROI %03d missing > .75 of voxels, exclude from analysis!!' %(r)
            command = """roi%03d=np.empty((timeslen,))""" %(r)
            exec(command)
            command = """roi%03d.fill(np.nan)""" %(r)
            exec(command)
        else: # Otherwise create timeseries within that ROI
            # Average timeseries within ROIs (masked) and output to a variable
            command = """3dmaskave -quiet -mask %s %s/%s""" %(roifile_masked,epidir,epi_filt)
            p=subprocess.Popen(command,stdout=subprocess.PIPE,shell=True)
            (tmpout,err) = p.communicate()
            out=tmpout.split()
            command = """roi%03d=np.asarray(out)""" %(r)
            exec(command)
            command = """roi%03d=roi%03d.astype(np.float)""" %(r,r) # Change from array of strings to floats
            exec(command)

        #add timeseries to full matrix
        if count==0: # Need to make the roi_all variable
            command = """roi_all=np.hstack((roi%03d[:,np.newaxis]))""" %(r)
        elif count==1: # Need to add dimension to roi_all b/c it's 1D
            command = """roi_all=np.hstack((roi_all[:,np.newaxis],roi%03d[:,np.newaxis]))""" %(r)
        else: # Add to existing roi_all variable
            command = """roi_all=np.hstack((roi_all,roi%03d[:,np.newaxis]))""" %(r)
        exec(command)
        
    roi_correlate=roi_all[start:stop,] # Include only the timepoints of the subblock
    
    # Save timeseries file across all ROIs
    timesoutfile='%s/%s' %(timesdir,timesoutname)
    np.savetxt(timesoutfile,roi_correlate,fmt='%.12f',delimiter='\t',newline='\n')

    # Correlate average timeseries across all ROIs and save correlation matrix output
    outcorr=np.corrcoef(roi_correlate.T) # Transpose so correlate ROIs instead of timepoints
    outcorr_masked=np.multiply(outcorr,mask) # Zero connections that are too short--this is really only for dosenbach_all
    corroutfile='%s/%s' %(roi2roidir,corroutname)
    np.savetxt(corroutfile,outcorr_masked,fmt='%.12f',delimiter='\t',newline='\n')
    
######### INPUTS #########

## User input
project = str(sys.argv[1]) # Note start at 1, because 0 is script name
atlas = str(sys.argv[2])
subblocks = str(sys.argv[3])
dosmooth = str(sys.argv[4])
dofilter = str(sys.argv[5])
doscrub = str(sys.argv[6]) 
docorr = str(sys.argv[7])

#examples:
#project = 'MegaRest.TMS' or 'despolarity'
#atlas = 'aal' or 'power'
#subblocks = #of subblocks to split each block into (must be even divisor)
#dosmooth = 'smooth' or 'nosmooth' or 'alreadysmooth'
#dofilter = 'filter' or 'nofilter' or 'alreadyfilter'
#doscrub = 'scrub' or 'noscrub' or 'both'
#docorr = 'corr' or 'nocorr'

# Parameters
fwhm=6
filt=[0.0009,0.08]
GSR=False # False is the default
maskedroi_prop=.25 # Minimum proportion of voxels in masked ROI as compared to unmasked ROI
mindist=20 # Minimum distance for connection in dosenbach_all atlas (load a mask that zeros out all connections shorter than mindist)


## Study-specific variables
homedir = "/home/despo/" + project + "/"
#scriptdir = "/home/despo/jrcohen/gitrepos/podNtools/" #CG: removed so we don't
#need dependence on user
if project == 'despolarity':
    datadir = homedir + "data/"
    tasks = ['rest']
    taskdirs = ['/rest/afni_rest_']
elif project == 'MegaRest.TMS':
    datadir = homedir + 'Data/'
    tasks = ['rest']
    taskdirs = ['/Analysis/Rest/afni_rest/']
nuisdir = homedir + "masks/nuisance/"
maskdir = homedir + "masks/native_" + atlas + "/"
maskdir_masked = homedir + "masks/native_" + atlas + "/masked/"

#CG - making names more generic below:
if subblocks == '1':
    analdir = homedir + "analyses/corr_Hallquist/"
elif subblocks == '2':
    analdir = homedir + "analyses/corr_2subblocks_Hallquist/"
elif subblocks == '3':
    analdir = homedir + "analyses/corr_3subblocks_Hallquist/"
elif subblocks == '6':
    analdir = homedir + "analyses/corr_6subblocks_Hallquist/"
elif subblocks == '12':
    analdir = homedir + "analyses/corr_12subblocks_Hallquist/"
    
if not os.path.exists(analdir):
    os.mkdir(analdir)
corrdir = analdir + "corr_" + atlas + "/"

if not os.path.exists(corrdir):
    os.mkdir(corrdir)
roi2roidir = corrdir + "/roi2roi/"
if not os.path.exists(roi2roidir):
    os.mkdir(roi2roidir)
timesdir = roi2roidir + "/timeseries/"
if not os.path.exists(timesdir):
    os.mkdir(timesdir)

if project == 'despolarity':
    subjects = ['201','202','203','204','205','206','207','208','209','210','211','212','213','214','215','216','217','218','219','220','221','222','223'] # All subjects
    #subjects = ['201','202','203','205','207','208','209','211','212','214','215','217','218','219','220','221','222','223'] # Follow-up only 
    #subjects = ['201'] # Short for trouble-shooting
elif project == 'MegaRest.TMS':
    #subjects = ['101_Sham'] #for testing
    #subjects =['101_Sham','102_Sham','103_Sham','104_Sham','105_Sham','106_Sham','107_Sham','108_Sham','109_Sham','110_Sham','111_Sham','112_Sham','113_Sham','114_Sham','116_Sham','118_Sham','119_Sham','120_Sham','121_Sham','122_Sham','123_Sham','124_Sham','125_Sham','126_Sham','127_Sham','128_Sham','129_Sham',
    #'101_Dlpfc','102_Dlpfc','103_Dlpfc','104_Dlpfc','105_Dlpfc','106_Dlpfc','107_Dlpfc','108_Dlpfc','109_Dlpfc','110_Dlpfc','111_Dlpfc','112_Dlpfc','113_Dlpfc','114_Dlpfc','116_Dlpfc','118_Dlpfc','119_Dlpfc','120_Dlpfc','121_Dlpfc','122_Dlpfc','123_Dlpfc','124_Dlpfc','125_Dlpfc','126_Dlpfc','127_Dlpfc','128_Dlpfc','129_Dlpfc',
    subjects = ['101_Aifo','102_Aifo','103_Aifo','104_Aifo','105_Aifo','106_Aifo','107_Aifo','108_Aifo','109_Aifo','110_Aifo','111_Aifo','112_Aifo','113_Aifo','114_Aifo','116_Aifo','118_Aifo','119_Aifo','120_Aifo','121_Aifo','122_Aifo','123_Aifo','124_Aifo','125_Aifo','126_Aifo','127_Aifo','128_Aifo','129_Aifo']


# Load parameter file for additional study-specific info
#params = np.loadtxt(scriptdir + '/funccon_params_' + project +
#'.txt',dtype=[('sub','|S20'),('nconds','int'),('nblocks','int'),('cond1','int'),('cond2','int'),('cond3','int')])
##removing dependence on script dir
params = np.loadtxt('funccon_params_' + project + '.txt',dtype=[('sub','|S20'),('nconds','int'),('nblocks','int'),('cond1','int'),('cond2','int'),('cond3','int')])

# Define atlases
if atlas == 'aal':
    roilist = range(1,91)
    #roilist = range(1,5) # Short for trouble-shooting
elif atlas == 'dosenbach':
    roilist = range(1,19) # Original 18 dosenbach ROIs
elif atlas == 'dosenbach_all' or atlas == 'power':
    roilist = range(1,265) # New dosenbach ROIs (Power paper)
elif atlas == 'HO' or atlas == 'HOall':
    roilist = range(1,97)
elif atlas == 'Greicius_func':
    roilist = range(1,91)
else:
    print 'This atlas has not been defined.'
    sys.exit()

    
######### MAIN SCRIPT #########
for sind,s in enumerate(subjects):

    print '      >>>> SUBJECT: %s' %s

    #establish the number (and order) of conditions for this subject
    if params['nconds'][np.where(params['sub']==s)]==3:
        conditions=[params['cond1'][np.where(params['sub']==s)],params['cond2'][np.where(params['sub']==s)],params['cond3'][np.where(params['sub']==s)]]
    elif params['nconds'][np.where(params['sub']==s)]==2:
        conditions=[params['cond1'][np.where(params['sub']==s)],params['cond2'][np.where(params['sub']==s)]]
    elif params['nconds'][np.where(params['sub']==s)]==1:
        conditions=[params['cond1'][np.where(params['sub']==s)]]
    elif params['nconds'][np.where(params['sub']==s)]==0:
        conditions=[]
    else:
        print 'This number of conditions has not been defined.'
        sys.exit()
    
    # Make conditions a list of strings
    for i in range(0,len(conditions)):
        conditions[i]=str(int(conditions[i]))
    if len(conditions) == 0: #no conditions
        conditions = [''] #empty string

    # Loop through conditions for each subject
    for cind,c in enumerate(conditions):

        #adjust directory name based on condition
        #CG: in the long run, it would probably be best to code this somewhere else?
        if project == 'MegaRest.TMS':
            taskdir_suff = ''
        else:
            taskdir_suff = s + c + '/' #added '/' here

        #loop through data (task) types
        for t in range(len(tasks)):
            epidir = datadir + s + taskdirs[t] + taskdir_suff
            bcount=1            

            #loop through blocks
            for b in range(1,params['nblocks'][np.where(params['sub']==s)]+1):

                # Define filenames
                # CG: these are somewhat project specific; specify in a separate place??
                if project == 'despolarity':
                    epi='%s%s-EPI-0%02d_proc.nii' %(s,c,b) 
                    epi_smooth='%s%s-EPI-0%02d_proc_smoothed.nii' %(s,c,b)
                    epi_filt='%s%s-EPI-0%02d_proc_filtered.nii' %(s,c,b)
                    regoutname='NUISANCE_ALL_r%02d.txt' %(b)
                elif project == 'MegaRest.TMS':
                    epi = 'pb02.%s.r%02d.volreg+orig' %(s,b) #motion corrected data
                    epi_smooth = 'pb03.%s.r%02d.smooth_Hallquist.nii' %(s,b) #smoothed data
                    epi_filt = 'pb04.%s.r%02d.filtered_Hallquist.nii' %(s,b) #filtered data
                    regoutname='NUISANCE_ALL_r%02d.txt' %(b)
                
                # First smooth epi data
                if dosmooth=='smooth':
                    smooth(epidir,epi,epi_smooth,fwhm)      
                
                # Next create nuisance regressor matrix
                if dofilter=='filter': # If nofilter or alreadyfilter don't need to get nuisance regressors or run bandpass filtering (means it was either already done, or only running correlations)
                    if dosmooth=='nosmooth':
                        get_nuisance_data(s,c,b,epidir,epi,nuisdir,regoutname,GSR)
                    else: # default is to use smoothed data
                        get_nuisance_data(s,c,b,epidir,epi_smooth,nuisdir,regoutname,GSR)
                        
                    # Now bandpass filter
                    if dosmooth=='nosmooth': 
                        do_bandpass_filtering(epidir,epi,epi_filt,regoutname,filt)
                    else: # default is to use smoothed data
                        do_bandpass_filtering(epidir,epi_smooth,epi_filt,regoutname,filt)
                    
                # Scrub data - CURRENTLY NOT CODED
                if doscrub == 'scrub' or doscrub == 'both':
                    scrub_data()
                
                # Last run regressions, separately for each run
                if docorr=='corr':
                    
                    # First determine length of timeseries (when dividing into
                    # N subblocks)
                    command = """fslnvols %s%s""" %(epidir,epi_filt)
                    p=subprocess.Popen(command,stdout=subprocess.PIPE,shell=True)
                    (tmpout,err) = p.communicate()
                    out=tmpout.split()
                    timeslen=int(out[0])
                    subblocklength=timeslen/int(subblocks)
                    
                    # Then define output names--using cind in case conditions aren't sequentially numbered for all subs
                    for subbl in range(1,int(subblocks)+1): 

                        #set up some counters
                        if subbl==1:
                            start=0
                            stop=subblocklength
                        else:
                            start=stop
                            stop=start+subblocklength

                        #filenames
                        if c == '':
                            corroutname='%s_Block%02d.txt' %(s,bcount) #bcount so cumulative across both blocks of a session
                            timesoutname='%s_Block%02d_timeseries.txt' %(s,bcount) #bcount so cumulative across both blocks of a session
                        else:
                            corroutname='%s_Sess%d_Block%02d.txt' %(s,cind+1,bcount) #bcount so cumulative across both blocks of a session
                            timesoutname='%s_Sess%d_Block%02d_timeseries.txt' %(s,cind+1,bcount) #bcount so cumulative across both blocks of a session
                        
                        
                        # Run correlations, selecting appropriate file
                        if dofilter=='nofilter':
                            if dosmooth=='nosmooth':
                                run_correlations(maskdir,maskdir_masked,atlas,roilist,maskedroi_prop,s,epidir,epi,timeslen,start,stop,roi2roidir,timesdir,corroutname,timesoutname,project) 
                            else: # Default is to use smoothed data
                                run_correlations(maskdir,maskdir_masked,atlas,roilist,maskedroi_prop,s,epidir,epi_smooth,timeslen,start,stop,roi2roidir,timesdir,corroutname,timesoutname,project)
                        else: # Default is to use filtered data
                            run_correlations(maskdir,maskdir_masked,atlas,roilist,maskedroi_prop,mindist,s,epidir,epi_filt,timeslen,start,stop,roi2roidir,timesdir,corroutname,timesoutname,project)
                            
                        ##### Have to incorporate options here for scrub (scrubbed data only) and both (both scrubbed and unscrubbed data)--for now just including all (unscrubbed data) #####

                        bcount=bcount+1
                        
