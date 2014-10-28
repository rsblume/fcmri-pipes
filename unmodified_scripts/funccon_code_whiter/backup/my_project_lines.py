# some project-specific code removed from merge:
# useful to setup config file!
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

# Get ROI names
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


