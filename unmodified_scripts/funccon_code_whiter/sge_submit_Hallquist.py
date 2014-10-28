#!/usr/bin/env python

# Stdlib
import os
from os.path import join as pjoin
import sys
import pandas

# Third-party
import numpy as np

command_file = '/home/despo/whiter/sge/sge_command_Hallquist.sh'
opt_file     = '/home/despo/whiter/sge/qsub_options.txt'

TEST = False
#TEST = True

#projects = ['bromo']
projects = ['pd']
#atlas = 'aal'
atlas = 'power'
suffix = 'rest'
subblock = '1'
# preproc options - don't do anything (debug)
#smooth_opt = 'nosmooth'
#filter_opt = 'nofilter'
#scrub_opt  = 'noscrub'
#corr_opt   = 'nocorr'

# preproc options - should make these subject-specific(?)
#smooth_opt = 'smooth'
#filter_opt = 'filter'
scrub_opt  = 'noscrub'
corr_opt   = 'corr'

# to fix missing ROIs or to run new atlas
smooth_opt = 'alreadysmooth'
filter_opt = 'alreadyfilter'

for project in projects:
# setup subject-specific params file
#from funccon_proc.py 
    basedir = '/home/despo/whiter/gitrepos/podNtools/funccon_code/params_files/'
    params_filename = basedir + 'funccon_params_' + suffix + '.' + project + '.csv' # Put name of this file into config file
    #print(params_filename)
    params_all = pandas.read_csv(params_filename,comment='#') # read in the subject-specific data
    params     = params_all.dropna(axis=0) # remove commented rows
    params[['cond','block']] = params[['cond','block']].astype(int) # final type int 
    for i in params.index:
        sub = params.subid[i]
        #this one is for funccon_bysubject.py
        command = 'qsub -@ %s -N fc.%s %s %s %s %s %s %s %s %s %s %s' % (opt_file, sub, command_file, project, atlas, suffix, subblock, smooth_opt, filter_opt, scrub_opt, corr_opt, i)
#ra!/usr/bin/env python
        print command
        if not TEST: # this is not a test
            os.system(command)

