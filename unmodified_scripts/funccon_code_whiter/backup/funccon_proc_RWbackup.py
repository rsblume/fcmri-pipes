#!/usr/bin/env python

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
import pandas

######### FUNCTIONS #########
def roilist_from_atlas(atlas):
    """ Returns a list of numeric ROI/node labels corresponding to the specificed atlas name
    """
    roilist_lookup_dict = { 
        'aal' : range(1,91),  # Automated anatomic labeling ROIs
         #'aal' : range(1,5) # Short for trouble-shooting
        'dosenbach' : range(1,19),  # Original 18 dosenbach ROIs
        'power'         : range(1,265), # Power paper (Neuron, 2011) ROIs
        'dosenbach_all' : range(1,265), # same as power 
        'HO'            : range(1,97),  # Harvard-Oxford ROIs
        'HOall'         : range(1,97),  # same as HO
        'Greicius_func' : range(1,91)  # Greicius paper (Shirer et al, Cereb Ctx. 2012?)
        }

    if roilist_lookup_dict.has_key(atlas):
       return roilist_lookup_dict[atlas]
    else:
        print 'The atlas ' + atlas + ' has not been defined.'
        sys.exit()

def ifisnii(epi):
    """ Function to determine if input epi file is an unzipped nifti file. This function assumes that this script must take .nii files as input. """

    if epi[-3:]=='nii':
        output=True
    else:
        print "Input EPI file must be .nii!!"
        sys.exit()
    
    return output


def define_epifilenames(epipath):
    """ Function to define filenames for input epi, smoothed epi, and filtered epi. This function assumes that input epi must be a .nii file. """

    filenames = dict()
    if ifisnii(params.epi[i]):
        filenames['epi']=params.epi[i].split('/')[-1] # This gets rid of the path so only has the epi name
        epidir=params.epi[i].rstrip(fi
        enames['epi']) # This removes epi chars from the end of the path only keep the path
        epibase=filenames['epi'][0:-4] # This gets rid of the suffix (.nii)
        filenames['epi_smooth']='%s_smoothed.nii' %(epibase)
        filenames['epi_filt']='%s_filtered.nii' %(epibase)

    return epidir,filenames

def remove_previous_file(filename):
    """ Function that forcibly removes previous versions of a file.
    """

    command = 'rm -f %s' %filename
    os.system(command)
    
def smooth(epidir,filenames,fwhm,dosmooth):
    ''' If dosmooth = smooth, smooth data and return filename of smoothed data.
    otherwise just return appropriate file name (smooth filename for
    alreadysmooth and original epi for nosmooth)

    Input:
    epidir: path (string) to directory where EPIs are stored
    filenames: dictionary of strings with EPI file naming conventions
    fwhm: value for smoothing (integer)
    dosmooth: string indicating whether to smooth or not ('smooth', 'nosmooth',
    or 'alreadysmooth')
    
    Return:
    filename (string) to use for the next step'''

    if dosmooth=='smooth':
        remove_previous_file('%s%s' %(epidir,filenames['epi_smooth'])) #added, otherwise doesn't resave
    
        #Note that in fox correlation script smoothing is equivalent to 3dmerge -1blur_fwhm; 3dBlurToFWHM better matches across subs so adopting it here
        command = """3dBlurToFWHM -input %s%s -prefix %s%s -FWHM %s""" %(epidir,filenames['epi'],epidir,filenames['epi_smooth'],fwhm)
        os.system(command)

        return filenames['epi_smooth']

    elif dosmooth == 'nosmooth':
        return filenames['epi']

    elif dosmooth == 'alreadysmooth':
        return filenames['epi_smooth']

    else:
        print 'Must define smooth parameter as smooth, nosmooth, or  alreadysmooth'
        sys.exit()

    
def get_nuisance_data(s,c,b,mot_gsr_epi_dir,nuis_dir,smooth_outfile,filenames,GSR=False):
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
            if project == 'rest.pd' or project == 'rest.bromo':
                # and thus, it begins... NEEDS FIXING
                roi_name='%s-EPI-%03d-1D.txt' % (s,b)
            else:
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
            print(command)
            os.system(command)

            # Derivative
            command = '1d_tool.py -infile %s/%s -derivative -demean -overwrite  -write %s/%s_deriv.r%02d.1D' %(epi_dir,roi_name,epi_dir,r,b)
            print(command)
            os.system(command)
            
        else: # Take the average intensity within the ROI first
            
            # Demean
            command = """3dmaskave -quiet -mask %s/%s %s/%s | 1d_tool.py -infile - -demean -overwrite -write %s/%s_demean.r%02d.1D""" %(roi_dir,roi_name,epi_dir,smooth_outfile,epi_dir,r,b)
            print(command)
            os.system(command)
            
            # Derivative
            command = """3dmaskave -quiet -mask %s/%s %s/%s | 1d_tool.py  -infile - -derivative -demean -overwrite -write %s/%s_deriv.r%02d.1D""" %(roi_dir,roi_name,epi_dir,smooth_outfile,epi_dir,r,b)
            print(command)
            os.system(command)
            
            
        # Now load individual regressor files so they can be concatenated
        command = """%s_demean=np.loadtxt('%s/%s_demean.r%02d.1D')""" %(r,epi_dir,r,b)
        print(command)
        exec(command)
        command = """%s_deriv=np.loadtxt('%s/%s_deriv.r%02d.1D')""" %(r,epi_dir,r,b)
        print(command)
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
    regoutfile='%s%s' %(epi_dir,filenames['regoutname'])
    np.savetxt(regoutfile,nuis_all,fmt='%.6f',delimiter='\t',newline='\n')

def do_bandpass_filtering(epidir,filenames,prev_step_file,filt,dofilter):
    ''' Do bandpass filtering, regression of nuisance variables, and despiking
    and output filtered file along with filename  if dofilter = 'filter'.
    Otherwise return appropriate outputfile name (filtered file for
    dofilter='alreadyfilter' and prev_step_file for dofilter='nofilter'
    
    Input:
    epidir: path (str) to EPI directory
    filenames: dictionary of strings with appropriate filenames for each
      processing step
    prev_step_file: string with the filename from the previous step (whatever
      is output as smooth_outfile currently) and will be operated on here/passed
      along if no filtering is supposed to happen.
    filt: array with low and high boundaries of bandpass filter
    dofilter: string specifying what to do ('filter', 'nofilter', or 'alreadyfilter'    
    '''

    if dofilter == 'filter':
        remove_previous_file('%s%s' %(epidir,filenames['epi_filt'])) #added, otherwise doesn't resave

        command = """3dBandpass -despike -ort %s%s -band %.4f %.4f -prefix %s%s -input %s%s""" %(epidir,filenames['regoutname'],filt[0],filt[1],epidir,filenames['epi_filt'],epidir,prev_step_file)
        os.system(command)

        return filenames['epi_filt']

    elif dofilter == 'nofilter':
        return prev_step_file

    elif dofilter == 'alreadyfilter':
        return filenames['epi_filt']

    else:
        print 'Must define filtering parameter as filter, nofilter, or alreadyfilter'
        sys.exit()
    
        
def scrub_data(doscrub,prev_step_file):
    ''' Scrub data to only include low-motion volumes '''
    ##### Holding off on this until determine the best way to execute scrubbing #####
    ##### Note that have to incorporate an option to do both scrubbed and
##### unscrubbed correlations in /"un_correlations function if doscrub=='both'
##### #####

    #Currently will simply return the appropriate filename
    if doscrub == 'scrub':
        print 'This section of the code still needs to be added!'
        raise NotImplementedError('scrub not implemented')
    elif doscrub == 'both':
        print 'This section of the code still needs to be added!'
        raise NotImplementedError('scrub and noscrub not implemented')
    elif doscrub == 'noscrub':
        #return the file from the previous step if no scrubbing happens.
        return prev_step_file
    else:
        print 'Must specify scrub parameter as scrub, both, or noscrub'
        raise ValueError('invalid scrub param: {0}'.format(doscrub))


def length_of_timeseries(epidir,filename,subblocks):
    """ Function for calculating the length of the time-series once divided
    into N subblocks
    """

    command = """fslnvols %s%s""" %(epidir,filename)
    p=subprocess.Popen(command,stdout=subprocess.PIPE,shell=True)
    (tmpout,err) = p.communicate()
    out=tmpout.split()
    timeslen=int(out[0])
    subblocklength=timeslen/int(subblocks)
    
    return subblocklength


        
def run_correlations(maskdir,maskdir_masked,atlas,roilist,maskedroi_prop,mindist,s,epidir,prev_output_file,timeslen,start,stop,roi2roidir,timesdir,filenames,project):
    ''' Run roi2roi correlations '''

    # First define and load anatomical distance mask if atlas is dosenbach_all
    # or power
    nnodes=len(roilist)
    if atlas=='dosenbach_all' or atlas == 'power':
        if project == 'despolarity' or project == 'task_graphtheory':
            maskfilename='%sanatomdist_masks/%s_%s_mindist%s.txt' %(roi2roidir,s,atlas,mindist)
        else:
            maskfilename = '%s%s_dist_over%s.txt' %(maskdir,atlas,mindist)
        
        if os.path.exists(maskfilename):
            mask=np.loadtxt(maskfilename)
        else:
            print '.... No distance file found... Assuming all connections can stay, regardless of distance'
            mask=np.ones((nnodes,nnodes)) # If doesn't exist, assume all connections can stay
    else:
        mask=np.ones((nnodes,nnodes)) # Also make a mask of 1s for other atlases to make multiplying easier below
        
    for count,r in enumerate(roilist):
        if atlas == 'aal' or atlas == 'shen':
            rstr='r%03d' %(r)
        elif atlas == 'Greicius_func' or atlas == 'dosenbach_all' or atlas == 'power':
            rstr='%03d' %(r)
        elif atlas == 'dosenbach' or atlas == 'HO':
            rstr='%02d' %(r)

        #get ROI file names. Clunky to have this if/else here. Should fix. 
        # JRC addition: Note that these names are identical except for the subject number (either the entire number or only the first 3 digits)
        # RW  note: using "+" to construct roinamei (line 318) fails if s/subid is passed as int!
        if project == 'MegaRest.TMS':
            roiname = '%s_%s_%s' %(atlas,rstr,s[:3])
            roifile = maskdir + roiname + 'r.nii'
            roifile_masked = maskdir_masked + roiname + '.nii'    
        elif project == 'rest.pd':
            roiname = '%s_%s_%sr.nii' %(atlas,rstr,s) # Note the trailing 'r' ...
            roifile = maskdir + roiname
            roiname = '%s_%s_%s.nii' %(atlas,rstr,s)
            roifile_masked = maskdir_masked + roiname
        elif project == 'rest.bromo':
            roiname = '%s_%s_%sr.nii.gz' %(atlas,rstr,s) # Note the trailing 'r' ...
            roifile = maskdir + roiname
            roiname = '%s_%s_%s.nii.gz' %(atlas,rstr,s)
            roifile_masked = maskdir_masked + roiname
        
        else:
            roiname =  atlas + "_" + rstr + "_" + s + "r.nii" 
            roifile = maskdir + roiname
            roifile_masked = maskdir_masked + roiname

        # Check to make sure enough of masked ROI exists before getting timeseries

        # Full ROI        
        # RW: Need to make sure these are available from prior preprocessing. How to fail gracefully if not?
        command = """fslstats %s -V""" %(roifile)
        #print(command)
        p=subprocess.Popen(command,stdout=subprocess.PIPE,shell=True)
        (tmpout,err) = p.communicate()
        out=tmpout.split()
        roivol=float(out[0]) # Make float so can divide and not output an integer
                
        # Masked ROI        
        command = """fslstats %s -V""" %(roifile_masked)
        #print(command)
        p=subprocess.Popen(command,stdout=subprocess.PIPE,shell=True)
        (tmpout,err) = p.communicate()
        out=tmpout.split()
        roivol_masked=float(out[0]) # Make float so can divide and not output an integer

        # Check proportion of voxels in masked ROI
        #if 0:
            #print("Warning: skipping check of ROI masking!")
        #print roivol_masked/roivol 
        if roivol_masked/roivol < maskedroi_prop:
            # Not enough of ROI left after masking; create a timeseries of nans
            print 'ROI %03d missing > .75 of voxels, exclude from analysis!!' %(r)
            command = """roi%03d=np.empty((timeslen,))""" %(r)
            exec(command)
            command = """roi%03d.fill(np.nan)""" %(r)
            exec(command)
        else: # Otherwise create timeseries within that ROI
            # Average timeseries within ROIs (masked) and output to a variable
            command = """3dmaskave -quiet -mask %s %s/%s""" %(roifile_masked,epidir,prev_output_file)
            print command
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
    timesoutfile='%s/%s' %(timesdir,filenames['timesoutname'])
    np.savetxt(timesoutfile,roi_correlate,fmt='%.12f',delimiter='\t',newline='\n')

    # Correlate average timeseries across all ROIs and save correlation matrix output
    outcorr=np.corrcoef(roi_correlate.T) # Transpose so correlate ROIs instead of timepoints
    outcorr_masked=np.multiply(outcorr,mask) # Zero connections that are too short--this is really only for dosenbach_all/power
    corroutfile='%s/%s' %(roi2roidir,filenames['corroutname'])
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
maskedroi_prop=.20 # Changed for rest.bromo 4061/4062 4/1/14
mindist=20 # Minimum distance for connection in dosenbach_all atlas (load a mask that zeros out all connections shorter than mindist)

## Study-specific variables
if project == 'despolarity':
    homedir = "/home/despo/" + project + "/"
    datadir = homedir + "data/"
    tasks = ['rest']
    taskdirs = ['/rest/afni_rest_']
    nuisdir = homedir + "masks/nuisance/"
    maskdir = homedir + "masks/native_" + atlas + "/"
    maskdir_masked = homedir + "masks/native_" + atlas + "/masked/"
elif project == 'MegaRest.TMS':
    homedir = "/home/despo/cgratton/data/" + project + "/"
    datadir = homedir + 'Data/'
    tasks = ['rest']
    taskdirs = ['/Analysis/Rest/afni_rest/']
    nuisdir = datadir + "Masks/nuisance/"
    maskdir = datadir + "Masks/native_" + atlas + "/"
    maskdir_masked = datadir + "Masks/native_" + atlas + "/masked/"
elif project == 'rest.pd':
    homedir = "/home/despo/pdwm/Rest/"
    datadir = homedir + 'Data/'
    nuisdir = datadir + "Masks/nuisance/masked/"
    maskdir = datadir + "Masks/native_" + atlas + "/"
    maskdir_masked = maskdir + "masked/"
elif project == 'rest.bromo':
    homedir = "/home/despo/rest.bromo/"
    datadir = homedir + 'Data/combined/'
    nuisdir = datadir + "Masks/nuisance/masked/"
    maskdir = datadir + "Masks/native_" + atlas + "/"
    maskdir_masked = maskdir + "masked/"
else:
    print 'MAKE SURE TO DEFINE datadir, tasks, taskdirs for this project!!!'
    sys.exit()


#CG - making names more generic below:
if project=='despolarity':
    if subblocks == '1':
        analdir = homedir + "/analyses/14_rest_corr_11-18-13_Hallquist/"
    elif subblocks == '2':
        analdir = homedir + "/analyses/16_rest_corr_2subblocks_11-18-13_Hallquist/"
    elif subblocks == '3':
        analdir = homedir + "/analyses/18_rest_corr_3subblocks_11-18-13_Hallquist/"
    elif subblocks == '6':
        analdir = homedir + "/analyses/20_rest_corr_6subblocks_11-18-13_Hallquist/"
    elif subblocks == '12':
        analdir = homedir + "/analyses/22_rest_corr_12subblocks_11-18-13_Hallquist/"
else:
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
    os.makedirs(analdir)
    # makedirs is like shell 'mkdir -p' : makes all intermediate-level directories needed to contain the leaf directory.
corrdir = analdir + "corr_" + atlas + "/"
if not os.path.exists(corrdir):
    os.mkdir(corrdir)
roi2roidir = corrdir + "/roi2roi/"
if not os.path.exists(roi2roidir):
    os.mkdir(roi2roidir)
timesdir = roi2roidir + "/timeseries/"
if not os.path.exists(timesdir):
    os.mkdir(timesdir)


# Load parameter file for additional study-specific info--must be in the same directory as this script!
#params = np.loadtxt('funccon_params_' + project + '.txt',dtype=[('sub','|S20'),('nconds','int'),('nblocks','int'),('cond1','int'),('cond2','int'),('cond3','int')])
params_filename = 'funccon_params_' + project + '.csv' # Put name of this file into config file
params_all = pandas.read_csv(params_filename,comment='#') # read in the subject-specific data
params     = params_all.dropna(axis=0) # remove commented rows
params[['subid','cond','block']] = params[['subid','cond','block']].astype(int) # final type int 

#### NOTE THAT NOW THIS HAS SUB, COND, BLOCK, EPI NAME. MAKE ANY VARIABLE 0 IF DON'T HAVE IT? ADD OTHERS PEOPLE MAY WANT? This has the dual benefit of having a readable list of all subjects (and relevant info) that have been run through this script, and not caring about specific study's session/condition/block/etc structure.


# Load atlas information
roilist = roilist_from_atlas(atlas)

######### MAIN SCRIPT #########

## Only loop needed is to loop through pandas file
#for i,idx in enumerate(params.subid): # Chose first column because they're all the same length--better way to access rows?
for i,irow in params.iterrows(): # pandas uses iterrows, now sure how best to use irow -- leave as params.<field>[i]
    print 'SUBID: ',params.subid[i],'\tCOND: ',params.cond[i],'\tBLOCK: ',params.block[i]
    epidir,filenames=define_epifilenames(params.epi[i])
    filenames['regoutname']='NUISANCE_ALL_r%02d.txt' %(params.block[i])
    ### NO LONGER NEED TO LOOP THROUGH SUBJECTS, CONDS, BLOCKS, NOR DEFINE THOSE VARIABLES OR SUBJECT-SPECIFIC DIRECTORIES

    # First smooth epi data
    smooth_outfile = smooth(epidir,filenames,fwhm,dosmooth)      
    #Note: could change smooth_outfile (and all following outfiles)
    # to prev_step_outfile to make the order of processing steps more generic in the future.
    #However, this seemed potentially more dangerous as well, so
    #not done for now.

    # Next create nuisance regressor matrix
    # If nofilter or alreadyfilter don't need to get nuisance regressors or run bandpass filtering (means it was either already done, or only running correlations)
    #OR that we're wanting to run correlations on unfiltered files? Maybe we should leave this an option to create either way?
    get_nuisance_data(str(params.subid[i]),str(params.cond[i]),
                      params.block[i],epidir,nuisdir,smooth_outfile,filenames,GSR)
                
    filter_outfile = do_bandpass_filtering(epidir,filenames,smooth_outfile,filt,dofilter)

    # Scrub data - CURRENTLY NOT CODED
    scrub_outfile = scrub_data(doscrub,filter_outfile)

    # Last run regressions, separately for each run
    if docorr=='corr':

        # First determine length of timeseries (when dividing into
        # N subblocks)
    	subblocklength = length_of_timeseries(epidir,scrub_outfile,subblocks)

        for subbl in range(1,int(subblocks)+1): 

            ##### NEED THIS IF GET RID OF COND AND BLOCK LOOPS AND USE JUST A FILE LIST!!!! #####
            # Define count variable for subblocks to name output files that resets for each new condition
            bcount=(params.block[i]-1)*int(subblocks)+subbl


            #set up some counters
            if subbl==1:
            	start=0
                stop=subblocklength
            else:
            	start=stop
                stop=start+subblocklength

            # Filenames
            if params.cond[i]==0: ## This assumes user uses "0" if they don't have a condition variable
                filenames['corroutname']='%s_Block%02d.txt' %(params.subid[i],bcount) #bcount so cumulative across both blocks of a session
                filenames['timesoutname']='%s_Block%02d_timeseries.txt' %(params.subid[i],bcount) #bcount so cumulative across both blocks of a session
            else:
                filenames['corroutname']='%s_Sess%d_Block%02d.txt' %(params.subid[i],params.cond[i],bcount) #bcount so cumulative across both blocks of a session
                filenames['timesoutname']='%s_Sess%d_Block%02d_timeseries.txt' %(params.subid[i],params.cond[i],bcount) #bcount so cumulative across both blocks of a session               

            # Run correlations, selecting appropriate file
            run_correlations(maskdir,maskdir_masked,atlas,roilist,maskedroi_prop,mindist,params.subid[i],epidir,scrub_outfile,subblocklength,start,stop,roi2roidir,timesdir,filenames,project) 
# REPLACED s WITH params.subid
            #note: changed timeslen to subblocklength here, since
            #that seemed like what you'd want? If that's incorrect,
            #change back.
                        
            ##### Have to incorporate options here for scrub (scrubbed data only) and both (both scrubbed and unscrubbed data)--for now just including all (unscrubbed data) #####

            #bcount=bcount+1 ## No longer needed here, define using params variable above
