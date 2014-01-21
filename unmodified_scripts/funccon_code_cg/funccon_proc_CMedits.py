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
import argparse
import tempfile

import numpy as np
import numpy.testing as npt

import pandas
import nibabel as ni

######### FUNCTIONS #########

def atlas_from_csv(infile):
    """ read atlas nodes from a csv file
    return as a list
    first line should be a header containing atlas name
    rest of the data should be in a columns"""
    dat = pandas.read_csv(infile, sep = None)
    return dat.values.aslist()


def atlas_init(atlas):
    """ return list of nodes based on the name of the atlas you are using,
    
    Parameters:
        atlas : string
            name of the atlas ('aal', doesnbach'...)
    Returns:
        roilist : list
            list of the nodes all numbering starts at 1
    """
    if atlas == 'aal':
        roilist = range(1,91)
        #roilist = range(1,5) # Short for trouble-shooting
    elif atlas == 'dosenbach':
        roilist = range(1,19) # Original 18 dosenbach ROIs
    elif atlas == 'dosenbach_all':
        roilist = range(1,265) # New dosenbach ROIs
    elif atlas == 'HO' or atlas == 'HOall':
        roilist = range(1,97)
    elif atlas == 'Greicius_func':
        roilist = range(1,91)
    else:
        raise ValueError('%s is not a valid atlas'%atlas)
    return roilist
   

def test_atlas_init():
    ### make sure atlas_initi gives me a reasonable return
    atlases = (
            'aal', 
            'dosenbach', 
            'HO', 
            'dosenbach_all', 
            'HOall',
            'Greicius_func')
    for atlas in atlases:
        mylist = atlas_init(atlas)
        ## all should return lists
        assert type(mylist) == type([])
        ## all should start with 1
        assert mylist[0] == 1 
    
    tmplist = atlas_init('aal')
    assert 90 in tmplist # 90 nodes in aal
    npt.assert_raises(ValueError, atlas_init, 'nonsense')

def smooth(epidir,epi,epi_smooth,fwhm):
    ''' Smooth EPI data '''

    # Note that in fox correlation script smoothing is equivalent to 3dmerge -1blur_fwhm; 3dBlurToFWHM better matches across subs so adopting it here
    command = """3dBlurToFWHM -input %s/%s -prefix %s/%s -FWHM %s""" %(epidir,epi,epidir,epi_smooth,fwhm)
    os.system(command)

def test_smooth():
    scriptdir = os.getcwd()

    tmpdir = tempfile.mkdtemp()
    tmpdat = np.zeros((10,10,10))
    tmpdat[5,5,5] = 1
    img = ni.Nifti1Image(tmpdat, np.eye(4))
    tmpfile = os.path.join(tmpdir, 'mytest.nii.gz')
    img.to_filename(tmpfile)

    ## move to dir to work on file
    os.chdir(tmpdir)
    ## run smooth
    os.chdir(startdir)


def get_nuisance_data(s,c,b,mot_gsr_epi_dir,epi_name,nuis_dir,regoutname,GSR=False):
    ''' Get mean and derivative of all motion parameters, white matter/CSF nuisance ROIs, and, optionally, whole brain signal '''

    if GSR==False:
        roi_type = ['mot','white','ventricles']
    else:
        roi_type = ['mot','white','ventricles','gsr']
    reg_type=['deriv','demean']

    for count,r in enumerate(roi_type):
        if r=='mot':
            epi_dir=mot_gsr_epi_dir
            roi_name='dfile.r%02d.1D' %(b)
        elif r=='white' or r=='ventricles':      
            roi_dir=nuis_dir
            roi_name='%s_%s.nii' %(r,s)
            epi_dir=mot_gsr_epi_dir
        elif r=='gsr':
            roi_dir=mot_gsr_epi_dir
            roi_name='%s-Mask-ALLSESSRUNS.nii' %(s)
            epi_dir=mot_gsr_epi_dir

        if r=='mot': # Just take the motion file output
            # Demean
            command = '1d_tool.py -infile %s/%s -demean -write %s/%s_demean.r%02d.1D' %(epi_dir,roi_name,epi_dir,r,b)
            os.system(command)
            # Derivative
            command = '1d_tool.py -infile %s/%s -derivative -demean -write %s/%s_deriv.r%02d.1D' %(epi_dir,roi_name,epi_dir,r,b)
            os.system(command)
        else: # Take the average intensity within the ROI first
            # Demean
            command = """3dmaskave -quiet -mask %s/%s %s/%s | 1d_tool.py -infile - -demean -write %s/%s_demean.r%02d.1D""" %(roi_dir,roi_name,epi_dir,epi_name,epi_dir,r,b)
            os.system(command)
            # Derivative
            command = """3dmaskave -quiet -mask %s/%s %s/%s | 1d_tool.py -infile - -derivative -demean -write %s/%s_deriv.r%02d.1D""" %(roi_dir,roi_name,epi_dir,epi_name,epi_dir,r,b)
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
    regoutfile='%s/%s' %(epi_dir,regoutname)
    np.savetxt(regoutfile,nuis_all,fmt='%.6f',delimiter='\t',newline='\n')

def do_bandpass_filtering(epidir,epi,epi_filt,regoutname,filt):
    ''' Do bandpass filtering, regression of nuisance variables, and despiking and output filtered file '''

    command = """3dBandpass -despike -ort %s/%s -band %.4f %.4f -prefix %s/%s -input %s/%s""" %(epidir,regoutname,filt[0],filt[1],epidir,epi_filt,epidir,epi)
    os.system(command)

def scrub_data():
    ''' Scrub data to only include low-motion volumes '''
    ##### Holding off on this until determine the best way to execute scrubbing #####
    ##### Note that have to incorporate an option to do both scrubbed and unscrubbed correlations in run_correlations function if doscrub=='both' #####

def run_correlations(maskdir,maskdir_masked,atlas,roilist,maskedroi_prop,mindist,s,epidir,epi_filt,timeslen,start,stop,roi2roidir,timesdir,corroutname,timesoutname):
    ''' Run roi2roi correlations '''

    # First define and load anatomical distance mask if atlas is dosenbach_all
    nnodes=len(roilist)
    if atlas=='dosenbach_all':
        maskfilename='%s/anatomdist_masks/%s_%s_mindist%s.txt' %(roi2roidir,s,atlas,mindist)
        if os.path.exists(maskfilename):
            mask=np.loadtxt(maskfilename)
        else:
            mask=np.ones((nnodes,nnodes)) # If doesn't exist, assume all connections can stay
    else:
        mask=np.ones((nnodes,nnodes)) # Also make a mask of 1s for other atlases to make multiplying easier below

    for count,r in enumerate(roilist):
        if atlas == 'aal':
            command="""rstr='r%03d'""" %(r)
            exec(command)
        elif atlas == 'Greicius_func' or atlas == 'dosenbach_all':
            command="""rstr='%03d'""" %(r)
            exec(command)
            rstr = '%03d'%r
        elif atlas == 'dosenbach' or atlas == 'HO':
            command="""rstr='%02d'""" %(r)
            exec(command)
        roiname =  atlas + "_" + rstr + "_" + s + "r.nii"
        roifile = maskdir + roiname
        roifile_masked = maskdir_masked + roiname


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


def run_main(project, atlas, subblocks, dosmooth, dofilter, doscru, docorr): 
    # Parameters
    ## Make into GLOBALS
    FWHM = 6
    FILT = [0.0009,0.08]
    GSR = False # False is the default
    MASKEDROI_PROP = .25 # minimum proportion of voxels in masked roi as compared to unmasked roi
    MINDIST = 20 # minimum distance for connection in dosenbach_all atlas (load a mask that zeros out all connections shorter than mindist)

    ## study-specific variables
    homedir = "/home/despo/" + project + "/"
    scriptdir = "/home/despo/jrcohen/gitrepos/podntools/"
    datadir = homedir + "/data/"
    if project == 'despolarity':
        tasks = ['rest']
        taskdirs = ['/rest/afni_rest_']
    nuisdir = homedir + "/masks/nuisance/"
    maskdir = homedir + "/masks/native_" + atlas + "/"
    maskdir_masked = homedir + "/masks/native_" + atlas + "/masked/"
    if subblocks == '1':
        analdir = homedir + "/analyses/14_rest_corr_11-18-13_hallquist/"
    elif subblocks == '2':
        analdir = homedir + "/analyses/16_rest_corr_2subblocks_11-18-13_hallquist/"
    elif subblocks == '3':
        analdir = homedir + "/analyses/18_rest_corr_3subblocks_11-18-13_hallquist/"
    elif subblocks == '6':
        analdir = homedir + "/analyses/20_rest_corr_6subblocks_11-18-13_hallquist/"
    elif subblocks == '12':
        analdir = homedir + "/analyses/22_rest_corr_12subblocks_11-18-13_hallquist/"

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

    subjects = ['201','202','203','204','205','206','207','208','209','210','211','212','213','214','215','216','217','218','219','220','221','222','223'] # all subjects
    #subjects = ['201','202','203','205','207','208','209','211','212','214','215','217','218','219','220','221','222','223'] # follow-up only 
    #subjects = ['201'] # short for trouble-shooting

    # load parameter file for additional study-specific info
    params = np.loadtxt(scriptdir + '/funccon_params_' + project + '.txt',dtype=[('sub','|s20'),('nconds','int'),('nblocks','int'),('cond1','int'),('cond2','int'),('cond3','int')])

    # define atlases
    ######### main script #########
    for sind,s in enumerate(subjects):
        if params['nconds'][np.where(params['sub']==s)]==3:
            conditions=[params['cond1'][np.where(params['sub']==s)],params['cond2'][np.where(params['sub']==s)],params['cond3'][np.where(params['sub']==s)]]
        elif params['nconds'][np.where(params['sub']==s)]==2:
            conditions=[params['cond1'][np.where(params['sub']==s)],params['cond2'][np.where(params['sub']==s)]]
        else:
            print 'this number of conditions has not been defined.'
            sys.exit()
        # make conditions a list of strings
        for i in range(0,len(conditions)):
            conditions[i]=str(int(conditions[i]))
        for cind,c in enumerate(conditions):
            taskdir_suff = s + c
            for t in range(len(tasks)):
                epidir = datadir + s + taskdirs[t] + taskdir_suff + "/"
                bcount=1
                for b in range(1,params['nblocks'][np.where(params['sub']==s)]+1):
                    # define filenames
                    epi='%s%s-epi-0%02d_proc.nii' %(s,c,b)
                    epi_smooth='%s%s-epi-0%02d_proc_smoothed.nii' %(s,c,b)
                    epi_filt='%s%s-epi-0%02d_proc_filtered.nii' %(s,c,b)
                    regoutname='nuisance_all_r%02d.txt' %(b)
                    # first smooth epi data
                    if dosmooth=='smooth':
                        smooth(epidir,epi,epi_smooth,fwhm)
                    # next create nuisance regressor matrix
                    if dofilter=='filter': # if nofilter or alreadyfilter don't need to get nuisance regressors or run bandpass filtering (means it was either already done, or only running correlations)
                        if dosmooth=='nosmooth':
                            get_nuisance_data(s,c,b,epidir,epi,nuisdir,regoutname,gsr)
                        else: # default is to use smoothed data
                            get_nuisance_data(s,c,b,epidir,epi_smooth,nuisdir,regoutname,gsr)
                        # now bandpass filter
                        if dosmooth=='nosmooth': 
                            do_bandpass_filtering(epidir,epi,epi_filt,regoutname,filt)
                        else: # default is to use smoothed data
                            do_bandpass_filtering(epidir,epi_smooth,epi_filt,regoutname,filt)
                    # scrub data
                    if doscrub == 'scrub' or doscrub == 'both':
                        scrub_data()
                    # last run regressions, separately for each run
                    if docorr=='corr':
                        # first determine length of timeseries
                        command = """fslnvols %s/%s""" %(epidir,epi_filt)
                        p=subprocess.popen(command,stdout=subprocess.pipe,shell=true)
                        (tmpout,err) = p.communicate()
                        out=tmpout.split()
                        timeslen=int(out[0])
                        subblocklength=timeslen/int(subblocks)
                        # then define output names--using cind in case conditions aren't sequentially numbered for all subs
                        for subbl in range(1,int(subblocks)+1): 
                            if subbl==1:
                                start=0
                                stop=subblocklength
                            else:
                                start=stop
                                stop=start+subblocklength
                            corroutname='%s_sess%d_block%02d.txt' %(s,cind+1,bcount) #bcount so cumulative across both blocks of a session
                            timesoutname='%s_sess%d_block%02d_timeseries.txt' %(s,cind+1,bcount) #bcount so cumulative across both blocks of a session
                        # run correlations
                            if dofilter=='nofilter':
                                if dosmooth=='nosmooth':
                                    run_correlations(maskdir,maskdir_masked,atlas,roilist,maskedroi_prop,s,epidir,epi,timeslen,start,stop,roi2roidir,timesdir,corroutname,timesoutname) 
                                else: # default is to use smoothed data
                                    run_correlations(maskdir,maskdir_masked,atlas,roilist,maskedroi_prop,s,epidir,epi_smooth,timeslen,start,stop,roi2roidir,timesdir,corroutname,timesoutname)
                            else: # default is to use filtered data
                                run_correlations(maskdir,maskdir_masked,atlas,roilist,maskedroi_prop,mindist,s,epidir,epi_filt,timeslen,start,stop,roi2roidir,timesdir,corroutname,timesoutname)
                            ##### have to incorporate options here for scrub (scrubbed data only) and both (both scrubbed and unscrubbed data)--for now just including all (unscrubbed data) #####
                            bcount=bcount+1

######### inputs #########
if __name__ == '__main__':
    ## user input
    parser = argparse.ArgumentParser()

    ## project
    parser.add_argument(
        'project',
        help = 'Name of the project'
        )
    parser.add_argument(
        'atlas',
        type = str, 
        help = 'atlas name, eg ("aal", "dosenbach")'
        )
    parser.add_argument(
        'subblocks',
        type = str, 
        help = """number of subblocks to split each block into--user's 
        responsibility to choose a number of subblocks that divides 
        evenly into the total block length"""
        )
    parser.add_argument(
        'dosmooth',
        type = str, 
        help = """one of : 'smooth' 'alreadysmooth' 'nosmooth'"""
        )
    parser.add_argument(
        'dofilter',
        type = str, 
        help = """ one of 'filter' 'alreadyfilter' 'nofilter'"""
        )
    parser.add_argument(
        'doscrub',
        type = str, 
        help = """one of : 'scrub' 'noscrub' 'both'"""
        )
    parser.add_argument(
        'docorr',
        type = str, 
        help = """ one of : 'corr' 'nocorr'"""
        )
    if len(sys.argv) < 2:
        parser.print_help()
    else:
        args = parser.parse_args()
        print args
