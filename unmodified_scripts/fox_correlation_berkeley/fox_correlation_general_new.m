function [] = fox_correlation_general_new(project,atlas,reps,TR,allsubj,win_roi,vox2vox,btwn_roi,partcorr,do_betacorr)

% Adapted for despolarity by JRC 5/4/12
% Adapted for kki_adhd by JRC 9/25/12
% Adapted for task_graphtheory by JRC 12/26/12
% Adapted for kki_autism by JRC 3/5/13
% Adapted for rest.nback by JRC 4/25/13
% despolarity usage: fox_correlation_general_new('despolarity','ATLAS',N,2,'n','n','n','n','n',0)
% kki_adhd usage:   fox_correlation_general_new('kki_adhd','ATLAS',1,2.5,'n','n','n','n','n',0)
% task_graphtheory usage: fox_correlation_general_new('task_graphtheory','ATLAS',N,2,'n','n','n','n','n',0)
% kki_autism usage:  fox_correlation_general_new('kki_autism','ATLAS',1,2.5,'n','n','n','n','n',0)
% rest.nback usage:  fox_correlation_general_new('rest.nback','dosenbach',1,2,'n','n','n','n','n',0)
%   project = main study dir below /home/despo
%   atlas = where ROIs came from--and what ROI dir is named
%   reps = how many blocks to split the rest run into (can be 1,2,3,6,12)
%   TR = TR rest data aquired with
%   allsubj = do all subjects at once? 'y' or 'n'; if say 'n' can then do 
%       multiple subjects at a time with multiple matlab windows
%   win_roi = do correlation of all voxels w/in each ROI? 'y' or 'n'
%   vox2vox = do correlations ROI to voxel or ROI to ROI; 'y' =
%       ROI to voxel; 'n' = ROI to ROI; 'b' = both
%    btwn_roi = do correlations from one voxel to other ROIs **NOT SURE 
%       ABOUT THIS, I SHOULD ASK!!!
%   partcorr = do partial correlation? 'y' or 'n'
%   do_betacorr = do beta correlations? 1 if y, 0 if no

%%% NOTE THAT I ONLY ADAPTED THIS SCRIPT FOR LOOPS I'M CURRENTLY USING, SO
%%% I MAY NEED TO ADAPT IT MORE IF I DO SOMETHING NEW!


% This m-file MAY produce data appropriate for full-brain correlation (not
% coherence) analyses.  The method is taken from
%   Fox et.al., PNAS, 7/5/2005
% The steps are as follows:
%
% 0) CURRENTLY THE WHOLE BRAIN MODE IS NOT NORMALIZED TO 1000
% 1) Load data.
% 2) Band-pass filter the time series (0.009 < f < 0.08)
% 3) Load nuisance parameters
%    a) rigid-body motion parameters
%    b) whole-brain mean signal, if desired
%    c) ventricular signal
%    d) white matter signal
% 4) Perform linear regression to remove the above and approximates of 
%       their temporal derivatives-- i.e. y = Bx, where y is the time 
%       series for a given voxel, x is the matrix of regressors, and 
%       E = y - Bx is the desired output (I believe).
% 5) Do correlation maps across the whole brain, using a seed region taken
%       from previous contrasts (e.g. R frontal ROI).
% 6) Determine the degrees of freedom = #time points/correction factor.
%       Correction factor from Jenkins & Watts, Spectral Analysis & Its
%       Applications, 1968 (Holden-Day, San Francisco).
% 7) Calculate z-scores, with std.dev = 1/sqrt(df - 3)
% 8) Compute performance correlations.
%
% Inputs: 1) fnames: the subject number (without -?? suffix).  If not
%                entered, the user will be prompted for it.
%
% Outputs: None to calling function

% FOR SOME REASON, DEFAULTS ARE NOT LOADED BY SPM_DEFAULT CALL.
% *************************************************
%%%%%%%%%%%%%%%%%%%%%%
addpath /home/despo/jrcohen/gitrepos/podNtools/ % To access some functions I need
scriptdir = '/home/despo/jrcohen/gitrepos/podNtools/';

%spm_defaults; % for older SPMs
spm_get_defaults;
p_smooth = 6;
f = [0.009 0.08];
maskthresh = .1;
mias_cf = 2.34;     % Was Mia's value of 6, but concerned it's changed.  Will
%   use Fox's value for now, but needs to be recomputed.
singleblock=0; %only make this true if computing lowTR for single blocks
disp(' ');

tstart=tic; % So can timestamp how long everything takes, with toc(tstart)
%%%%%%%%%%%%%%%%%%%%

if isequal(project,'kki_adhd'),
    proc = input('Is this data already processed (1 if yes, 0 if no)? ');
    task = input('Is this data RS only or does it include GNG (1 if everything, 0 if RS only)? ');
end;

if isequal(project,'kki_autism'),
    proc = input('Is this data already processed (1 if yes, 0 if no)? ');
end;

ScriptPath = ['/home/despo/jrcohen/gitrepos/podNtools/'];
if isequal(project,'kki_adhd') && proc,
    StudyPath = ['/home/despo/' project '/data_processed/'];
elseif isequal(project,'kki_autism'),
    if proc,
        StudyPath = ['/home/despo/jrcohen/' project '/data_processed/'];
    else,
        StudyPath = ['/home/despo/jrcohen/' project '/data/'];
    end;
else,
    StudyPath = ['/home/despo/' project '/data/'];
end;
NuisPath_btwn = ['/home/despo/' project '/masks/nuisance_btwn/masked/']; % This doesn't exist for despolarity as of 5/4/12
NuisPath = ['/home/despo/' project '/masks/nuisance/']; % Note the original line had a /masked/ subfolder, but I din't mask my nuisance ROIs so I don't have that as of 5/4/12
% Add loop defining name of directory based on number of reps
if isequal(project,'despolarity'),
    if reps == 1,
        AnalPath = ['/home/despo/' project '/analyses/02_rest_corr_8-12-12/'];
    elseif reps == 2,
        AnalPath = ['/home/despo/' project '/analyses/04_rest_corr_2subblocks_9-15-12/'];
    elseif reps == 3,
        AnalPath = ['/home/despo/' project '/analyses/06_rest_corr_3subblocks_9-15-12/'];
    elseif reps == 6,
        AnalPath = ['/home/despo/' project '/analyses/08_rest_corr_6subblocks_9-15-12/'];
    elseif reps == 12,
        AnalPath = ['/home/despo/' project '/analyses/10_rest_corr_12subblocks_9-15-12/'];
    end;
elseif isequal(project,'kki_adhd'),
    if proc,
        if task && reps==1,
            AnalPath = ['/home/despo/' project '/analyses/01_rest_corr_2-5-13_KKIproc/'];
        elseif ~task && reps==1,
            AnalPath = ['/home/despo/' project '/analyses/03_rest_corr_4-22-13_KKIproc_restonly/'];
        else,
            fprintf('NEED TO DEFINE ANAL DIR!\n');
            return;
        end;
    else,
        fprintf('NEED TO DEFINE ANAL DIR!\n');
        return;
    end;
elseif isequal(project,'task_graphtheory'),
    if reps == 1,
        AnalPath = ['/home/despo/' project '/analyses/02_rest_corr_12-26-12/'];
    elseif reps == 2,
        fprintf('NEED TO DEFINE ANAL DIR!\n');
        return;
    elseif reps == 3,
        fprintf('NEED TO DEFINE ANAL DIR!\n');
        return;
    elseif reps == 6,
        fprintf('NEED TO DEFINE ANAL DIR!\n');
        return;
    elseif reps == 12,
        fprintf('NEED TO DEFINE ANAL DIR!\n');
        return;
    end;
elseif isequal(project,'kki_autism'),
    if proc,
        if reps==1,
            AnalPath = ['/home/despo/jrcohen/' project '/analyses/01_rest_corr_3-5-13_KKIproc/'];
        else,
            fprintf('NEED TO DEFINE ANAL DIR!\n');
            return;
        end;
    else,
        fprintf('NEED TO DEFINE ANAL DIR!\n');
        return;
    end;
elseif isequal(project,'rest.nback'),
    if reps==1,
        AnalPath = ['/home/despo/jrcohen/' project '/analyses/11_timeseries_corr_4-10-13/'];
    else,
        fprintf('NEED TO DEFINE ANAL DIR!\n');
        return;
    end;
else,
    fprintf('NEED TO DEFINE ANAL DIR!\n');
    return;
end;

if isequal(project,'Rest.TMS') || (isequal(project,'Rest.ShamTMS') && isequal(atlas,'aalhalf'))
    ROIPath = [StudyPath 'Masks/native_' atlas '/'];
elseif isequal(project,'despolarity') || isequal(project,'task_graphtheory')  || isequal(project,'rest.nback')
    ROIPath = ['/home/despo/' project '/masks/native_' atlas '/'];
elseif isequal(project,'kki_adhd'),
    if proc,
        ROIPath = ['/home/despo/' project '/masks/norm_' atlas '/'];
    else,
        ROIPath = ['/home/despo/' project '/masks/native_' atlas '/'];
    end;
elseif isequal(project,'kki_autism'),
    if proc,
        ROIPath = ['/home/despo/jrcohen/' project '/masks/norm_' atlas '/'];
    else,
        ROIPath = ['/home/despo/jrcohen/' project '/masks/native_' atlas '/'];
    end;
elseif isequal(atlas,'face_scene') || isequal(atlas,'what_where')
    ROIPath = [StudyPath 'Masks/native_' atlas '/epispace/'];
elseif isequal(project,'MegaRest.TMS_ASLconn')
    ROIPath = [StudyPath 'Masks/native_' atlas '/'];
elseif isequal(btwn_roi,'y')
    ROIPath = [StudyPath 'Masks/native_' atlas '/singlevox/'];
else
    ROIPath = [StudyPath 'Masks/native_' atlas '/masked/'];
end
if isequal(project,'despolarity') || isequal(project,'kki_adhd') || isequal(project,'task_graphtheory')  || isequal(project,'rest.nback')
    ROIPath_btwn = ['/home/despo/' project '/masks/native_' atlas '_singlevox/']; % This doesn't exist for JRC projects
else
    ROIPath_btwn = [StudyPath 'Masks/native_' atlas '_singlevox/'];
end

%LesionDir = [StudyPath 'Masks/native_lesion/lesion_epispace_singlevox/']; % Unecessary for JRC projects

output = input('What type of output do you want z/r/b (for graph theory want r): ','s'); 

if isequal(win_roi,'y')
    StudyPathnm = ['corr_' atlas '/'];
    smooth = 'n';
    disp('Not smoothing');
    ROIPath_suff = 'peel_';
elseif isequal(btwn_roi,'y')
    StudyPathnm = ['corr_' atlas '/'];
    smooth = 'n';
    disp('Not smoothing');
    ROIPath_suff = '';
    btwn_roi_val = input('Btwn which ROI?  ','s');
elseif isequal(output,'r')
    StudyPathnm = ['corr_' atlas '/'];
    if (isequal(project,'kki_adhd') || isequal(project,'kki_autism')) && proc,
        smooth = 'n';
    else,
        smooth = 'y';
        disp('smoothing');
    end;
    ROIPath_suff = '';
elseif isequal(output,'b')
    StudyPathnm_z = ['corr_z_' atlas '/'];
    StudyPathnm_r = ['corr_' atlas '/'];
    if isequal(project,'kki_adhd') && proc,
        smooth = 'n';
    else,
        smooth = 'y';
        disp('smoothing');
    end
    ROIPath_suff = '';
    StudyPathnm = StudyPathnm_z;
else
    StudyPathnm = ['corr_z_' atlas '/'];
    if isequal(project,'kki_adhd') && proc,
        smooth = 'n';
    else,
        smooth = 'y';
        disp('smoothing');
    end
    ROIPath_suff = '';
end
if do_betacorr
    disp('running beta correlations');
    StudyPathnm = ['corr_beta/']%['corr_beta_' atlas '/'];
    beta_cond = input('Which condition? (face_right_cor): ', 's');
end

EPIDIR = StudyPath;

dodosen = 0;
sessions = {'1'}; % Default, change in individual project loop below if necessary
if isequal(project,'Rest.ShamTMS')
    d = [dir([EPIDIR '/2*']);dir([EPIDIR '/1*'])];
    EPIDIR_suff = '/Rest/'; %NEW analysis (old used to be Total/NIfTI
    motionfilename = '*00*-1D.txt';
    epi_name_1 = 'EPI-0*';
    epi_name_2 = '*-CoReg.nii';
    regress_nuis = 1;
    blocklist = '123456';
    dodosen = 1;
elseif isequal(project,'Rest.OlderControl') || isequal(project,'Rest.MidControl')
    d = [dir([EPIDIR '/2*']);dir([EPIDIR '/1*'])];
    EPIDIR_suff = '/Total/NIfTI/';
    motionfilename = '*00*-1D.txt';
    epi_name_1 = 'EPI-0*';
    epi_name_2 = '*-CoReg.nii';
    regress_nuis = 1;
    blocklist = '123456';
    dodosen = 1;
elseif isequal(project,'Rest.thetaTMS')
    d = dir([EPIDIR '/3*']);
    prepost = input('Pre or Post? (pre/post): ', 's');
    EPIDIR_suff = ['/' prepost 'TMS_Rest/Analysis/afni_rest/'];
    motionfilename = 'dfile.r*.1D';
    epi_name_1 = 'EPI-0*';
    epi_name_2 = '*.nii';
    regress_nuis = 1;
    if isequal(prepost,'pre')
        blocklist = '1234';
    else
        blocklist = '12345678';
    end
elseif isequal(project,'Rest.TMS')
    d = dir([EPIDIR '/1*']);
    EPIDIR_suff = '/Total/NIfTI/';
    motionfilename = '*00*-1D.txt';
    epi_name_1 = 'EPI-0*';
    epi_name_2 = '*-CoReg.nii';
    regress_nuis = 1;
    blocklist = '123456';
elseif isequal(project,'Rest.Lesion')
    d = [dir([EPIDIR '/1*']); dir([EPIDIR '/3*'])]; %all the patients
    %EPIDIR_suff = '/Total/NIfTI/';
    EPIDIR_suff = '/Rest/';
    motionfilename = '*00*-1D.txt';
    %epi_name_1 = 'EPI-00*';
    epi_name_1 = 'EPI-0*';
    if isequal(btwn_roi,'y')
        epi_name_2 = '*-CoReg_resam.nii';
    else
        epi_name_2 = '*-CoReg.nii';
    end
    regress_nuis = 1;
    blocklist = '12';
elseif isequal(project,'Rest.Cerebellum')
    d = [dir([EPIDIR '/3*']); dir([EPIDIR '/4*']); dir([EPIDIR '/6*'])];
    EPIDIR_suff = '/Total/NIfTI/';
    motionfilename = '*00*-1D.txt';
    epi_name_1 = 'EPI-00*';
    epi_name_2 = '*-CoReg.nii';
    regress_nuis = 1;
    blocklist = '123456';
elseif isequal(project,'Atten.TMS') || isequal(project,'Atten.TMS2')
    d = dir([EPIDIR '/3*']);
    EPIDIR_suff = '/Lateralize/Analysis/afni_proc_block/';
    motionfilename = 'dfile.r*.1D';
    if do_betacorr
        epi_name_1 = ['beta.' beta_cond '.*'];%'beta.*_cor.*';
        epi_name_2 = '*.nii';
        regress_nuis = 0;
    else
        epi_name_1 = 'EPI-00*';
        epi_name_2 = '*.nii';
        regress_nuis = 1;
        regress_motion = 1
    end
    blocklist = '123456';
elseif isequal(project,'Rest.Norway')
    d = [dir([EPIDIR '/0*']); dir([EPIDIR '/1*'])];
    EPIDIR_suff = '/Analysis/Rest/afni_rest';
    motionfilename = 'dfile.r*.1D';
    epi_name_1 = 'EPI-0*';
    epi_name_2 = '*.nii';
    regress_nuis = 1;
    blocklist = '1';
elseif isequal(project,'Rest.LowTR')
    d = [dir([EPIDIR '/0*']); dir([EPIDIR '/1*'])];
    EPIDIR_suff = '/Analysis/Rest/afni_rest';
    motionfilename = 'dfile.r*.1D';
    epi_name_1 = 'EPI-00*';
    epi_name_2 = '*.nii';
    regress_nuis = 1;
    blocklist = '123456';
elseif isequal(project,'Rest.TMSnew')
    d = dir([EPIDIR '/1*']);
    EPIDIR_suff = '/Analysis/Rest/afni_rest';
    motionfilename = 'dfile.r*.1D';
    epi_name_1 = 'EPI-0*';
    epi_name_2 = '*.nii';
    regress_nuis = 1;
    blocklist = [1:10];%'123456';
elseif isequal(project,'MegaRest.TMS_ASLconn')
    d = dir([EPIDIR '/1*']);
    %EPIDIR_suff = '/Analysis/Rest/afni_rest';
    EPIDIR_suff = '/Analysis/ASL/';
    %motionfilename = 'dfile.r*.1D';
    motionfilename = 'motion_params*.txt';
    %epi_name_1 = 'EPI-0*';
    epi_name_1 = 'CBF_';
    %epi_name_2 = '*.nii';
    epi_name_2 = '*_tseries.nii';
    regress_nuis = 1;
    blocklist = [1:3];%'123456';
    NuisPath = [StudyPath 'Masks/nuisance/'];
    f = [0.009 0.06]; %slightly lower upper bound on freq due to low TR 
elseif isequal(project,'MegaRest.TMS')
    d = dir([EPIDIR '/1*']);
    EPIDIR_suff = '/Analysis/Rest/afni_rest';
    motionfilename = 'dfile.r*.1D';
    epi_name_1 = 'EPI-0*';
    epi_name_2 = '*.nii';
    regress_nuis = 1;
    blocklist = [1:3];%'123456';
    NuisPath = [StudyPath 'Masks/nuisance/'];
    %dodosen = 1; %CG - not sure that we needed this??
elseif isequal(project,'despolarity')
    d = dir([EPIDIR '/2*']);
    EPIDIR_suff = '/rest/afni_rest_';
    %sessions={'1','2','3'}; % Define this below because not all subjects had the same sessions
    motionfilename = 'dfile.r*.1D';
    epi_name_1 = 'EPI-00*';
    epi_name_2 = '*_proc.nii';
    regress_nuis = 1;
    blocklist = '12'; 
    NuisPath = ['/home/despo/' project '/masks/nuisance/'];
    dodosen = 0;
elseif isequal(project,'kki_adhd')
    tasks={'RS','GNG','GNG'};
    d = dir([EPIDIR '/' tasks{1} '/*']); % Take sub list from rest dir b/c that includes all possible subs
    if proc,
        EPIDIR_suff = {'/','/Functionals/','/Functionals/'};
        epi_name_1 = '';
        epi_name_2 = '_KKIproc.nii';
        tasktype = {'RS','GNGs','GNGr'};
        regress_mot=0;
        regress_nuis=0;
    else,
        EPIDIR_suff = '/rest/afni_rest_';
        motionfilename = 'dfile.r*.1D';
        epi_name_1 = 'EPI-00*';
        epi_name_2 = '*_proc.nii';
        regress_nuis = 1;
    end;
    sessions = {'1'};
    if task,
        blocklist = '123'; % RS,GNGs,GNGr 
    else,
        blocklist = '1'; 
    end;
    NuisPath = ['/home/despo/' project '/masks/nuisance/'];
    dodosen = 0;
elseif isequal(project,'task_graphtheory')
    d = dir([EPIDIR '/3*']);
    task = input('Which task? (rest, nback, bimancoord): ','s');
    while ~isequal(task,'rest') && ~isequal(task,'nback') && ~isequal(task,'bimancoord'),
        task = input('You must enter: rest, nback, bimancoord. Please re-enter: ','s');
    end;
    eval(sprintf('EPIDIR_suff = ''/%s/afni_%s_'';',task,task));
    if isequal(task,'rest'),
        concat=0; % No loads for rest so never concatenate
        sessions={'1','1','2','2'};
        normtype={'A','B','A','B'};
        epi_name_1 = 'EPI-00*';
        epi_name_2 = '*_proc.nii';
        blocklist = '12';
        blocksuff = {'',''};
    else
        concat = input('Are the task runs by block or concatenated by load? 1 if concat by load, 0 if using blocks.\n(note that it is better to concatenate AFTER, which is option 0): ');
        while isempty(find(concat==[0 1])),
            concat = input('You must enter 1 if concatenated by load, 0 if using blocks: ');
        end;
        %sessions={}; % Actual session with 1 or 2 defined below in main sub loop
        normtype={'A','B'};
        epi_name_2 = '*_proc.nii';
        if concat,
            epi_name_1 = 'EPI-';
            if isequal(task,'nback'),
                blocklist = '023';
                blocksuff = {'back','back','back'};
            elseif isequal(task,'bimancoord'),
                blocklist = 'ri';
                blocksuff = {'l','nt'};
            end
        else,
            epi_name_1 = 'EPI-00*';
            blocklist  = '1234';
            blocksuff = {'','','',''};
        end;
    end
    maskdir_pref = '/rest/afni_rest_';
    maskdir_suff='1';
    motionfilename = 'dfile.r*';
    motionfilename_suff = '.1D';
    regress_nuis = 1;
    NuisPath = ['/home/despo/' project '/masks/nuisance/'];
    dodosen = 0;
elseif isequal(project,'kki_autism')
    d = dir([EPIDIR '/*']);
    if proc,
        EPIDIR_suff = '/rest/';
        epi_name_1 = '';
        epi_name_2 = '.fsnwcfwa.nii';
        regress_mot=0;
        regress_nuis=0;
    else,
        EPIDIR_suff = '/rest/afni_rest_';
        motionfilename = 'dfile.r*.1D';
        epi_name_1 = 'EPI-00*';
        epi_name_2 = '*_proc.nii';
        regress_nuis = 1;
    end;
    sessions = {'1'};
    blocklist = '1'; 
    NuisPath = ['/home/despo/jrcohen/' project '/masks/nuisance/'];
    dodosen = 0;
elseif isequal(project,'rest.nback')
    d = dir([EPIDIR '/1*']);
    task = input('Which task? (rest, Nback, TRSE): ','s');
    while ~isequal(task,'rest') && ~isequal(task,'Nback') && ~isequal(task,'TRSE'),
        task = input('You must enter: rest, Nback, or TRSE. Please re-enter: ','s');
    end;
    if isequal(task,'rest'),
        taskdir='/rest/afni_';
        taskdirsuff='/';
    else,
        taskdir='/functional/afni_func_';
        taskdirsuff='_timeseries/';
    end;
    eval(sprintf('EPIDIR_suff = ''/%s%s%s'';',taskdir,task,taskdirsuff));
    epi_name_1 = 'EPI-0';
    epi_name_2 = '_proc.nii';
    if isequal(task,'rest'),
        blocklist = [1]; % Needs to be double because TRSE has > 10
        TR=1; % Define TR here because different for different tasks
    elseif isequal(task,'Nback'),
        blocklist = [1:2]; % Needs to be double because TRSE has > 10
        TR=2; % Define TR here because different for different tasks
    elseif isequal(task,'TRSE'),
        blocklist = [1:16]; % Needs to be double because TRSE has > 10
        TR=1; % Define TR here because different for different tasks
    end
    motionfilename = 'dfile.*';
    motionfilename_suff = '.1D';
    regress_nuis = 1;
    NuisPath = ['/home/despo/' project '/masks/nuisance/'];
    dodosen = 0;
end


wb_or_n = 'n';
if wb_or_n == 'y'
    StudyPathnm = ['corr_' atlas '_globalsignal/'];
    disp('HAVE A WHOLE-BRAIN REGRESSOR!!');
else
    disp(' NO whole-brain regressor.');
end


cf_or_n = 'c';
disp('  Artificially setting use Mia''s correction factor of ~six.');
while cf_or_n ~= 's' && cf_or_n ~= 'c'
    cf_or_n = input('Use Mia''s correction factor of ~six, or calculate (s/c) ?  ','s');
end

block_select = input('Do all blocks? (y/n): ', 's');
if isequal(block_select,'n')
    blocklist = input('Which block(s)? "1":  ','s');
end

if isequal(vox2vox,'y') || isequal(vox2vox,'b')
    
    mask_epi_data = input('Mask EPI data?: y/n ', 's');
    
    if isequal(mask_epi_data,'y') % Not tested for JRC projects
        if isequal(project,'despolarity') || isequal(project,'kki_adhd') || isequal(project,'task_graphtheory')  || isequal(project,'rest.nback')
            MaskPath = ['/home/despo/' project '/masks/native_' atlas '_PFC_mask/'];
        else
            MaskPath = [StudyPath 'Masks/native_' atlas '_PFC_mask/'];
        end
        mask_epi_data_file = ['lPFC';'mPFC'];
        hemi_type = input('Which mask hemispheres? L/R/B: ','s');
        if isequal(hemi_type,'B')
            hemi_list = ['L';'R'];
        else
            hemi_list = hemi_type;
        end
    end
end

done = 'n';

count=0;

if allsubj~='y'
    while done ~= 'y'
        fnames = input('Please enter the subject number(s) (eg: {''301A''} or {''301A'',''302A''}): ');
        for file=1:length(fnames)
            for l=1:length(d)
                if strcmp(d(l).name, fnames{file}) % strcmp instead of strfind so I have to type the exact names
                    count=count+1;
                    named(count).name=d(l).name;
                    named(count).fnum=fnames{file};
                    
                    if isequal(d(l).name,'01') ||isequal(d(l).name,'02') || isequal(d(l).name,'03') || isequal(d(l).name,'04') || isequal(d(l).name,'11')
                        singleblock = 1;
                    end
                end
            end
        end
        done = input('     Is the number correct (y/n) ?  ','s');
    end
else
    for l=1:length(d)
        count=count+1;
        named(count).name=d(l).name;
        named(count).fnum=d(l).name; %is this right?
    end
end


for l=1:length(named)
    
    if isequal(atlas,'face_scene') || isequal(atlas,'what_where')
        r=dir([ROIPath, ROIPath_suff '*_',named(l).name '_r.nii']);
        subname = named(l).name;
    elseif isequal(project,'MegaRest.TMS')
        r=dir([ROIPath, ROIPath_suff '*_*',named(l).name(1:3) '.nii']);
        subname = named(l).name(1:3);
    elseif isequal(project,'MegaRest.TMS_ASLconn')
        r=dir([ROIPath, ROIPath_suff 'r*_*',named(l).name(1:3) '.nii']);
        subname = named(l).name(1:3);
    elseif (isequal(project,'kki_adhd') || isequal(project,'kki_autism')) && proc,
        r=dir([ROIPath, ROIPath_suff '*r.nii']);
        subname = named(l).name;
    elseif isequal(project,'rest.nback') 
        r=dir([ROIPath, ROIPath_suff '*_*_' named(l).name 'r_' task '.nii']); 
        subname = named(l).name;
    else
        r=dir([ROIPath, ROIPath_suff '*_*_' named(l).name 'r.nii']); 
        subname = named(l).name;
    end
    
    if isequal(project,'despolarity'),
        if isequal(subname,'217'), 
            sessions={'1','3','4'}; % This subject had 4 sessions, but not using 2
        elseif isequal(subname,'215'),  
            sessions={'1','2','4'}; % This subject had 4 sessions, but not using 3   
        elseif isequal(subname,'204') || isequal(subname,'206') || isequal(subname,'210') || isequal(subname,'213') || isequal(subname,'216'),
            sessions={'1','2'}; % These subjects had only 2 sessions (no follow-up)
        else,
            sessions={'1','2','3'}; % The rest of the subjects had 3 sessions (2 intitial plus 1 follow-up)
        end;
    end;

    if isequal(atlas,'HOall') && length(r) < 96
        dbstop;
    end
    
    count=1;
    for l2=1:length(r)
        if ~isempty(strfind(r(l2).name, subname)) || ((isequal(project,'kki_adhd') || isequal(project,'kki_autism')) && proc)
            thisv = spm_vol([ROIPath, r(l2).name]);
            [thisY, thisXYZ] = spm_read_vols(thisv);
            named(l).roi(count).name = r(l2).name;
            named(l).roi(count).idx = find(thisY);
            count = count + 1;
            clear thisY thisXYZ;
        else
            disp('did not load any rois');
        end
    end
end

if isequal(btwn_roi,'y')
    for l=1:length(named)
        r = dir([ROIPath_btwn, ROIPath_suff '*' atlas btwn_roi_val '_',named(l).name '*.nii']);
        count=1;
        for l2=1:length(r)
            if strfind(r(l2).name, named(l).name)
                thisv = spm_vol([ROIPath_btwn, r(l2).name]);
                [thisY, thisXYZ] = spm_read_vols(thisv);
                named(l).btwnroi(count).name = r(l2).name;
                named(l).btwnroi(count).idx = find(thisY);
                count = count + 1;
                clear thisY thisXYZ;
            else
                disp('did not load any rois');
            end
        end
    end
    if regress_nuis
        for l=1:length(named)
            r=dir([NuisPath_btwn, '*',named(l).name '*.nii']);
            ncount=1;
            for l2=1:length(r)
                if ~isempty(strfind(r(l2).name, 'white')) | ~isempty(strfind(r(l2).name, 'ventricle'))
                    thisv = spm_vol([NuisPath_btwn, r(l2).name]);
                    [thisY, thisXYZ] = spm_read_vols(thisv);
                    named(l).nuisroi(ncount).name = r(l2).name;
                    named(l).nuisroi(ncount).idx = find(thisY);
                    ncount = ncount + 1;
                else
                    disp('did not load any nuisance rois');
                end
            end
        end
    end
else
    if regress_nuis
        for l=1:length(named)
            if isequal(project,'MegaRest.TMS')
                r=dir([NuisPath, '*',named(l).name(1:3) '*.nii']);
            elseif isequal(project,'MegaRest.TMS_ASLconn')
                r=dir([NuisPath, 'r*',named(l).name(1:3) '*.nii']);
            elseif isequal(project,'rest.nback')
                r=dir([NuisPath, '*',named(l).name,'*',task,'*.nii']);
            else
                r=dir([NuisPath, '*',named(l).name,'*.nii']);
            end
            ncount=1;
            for l2=1:length(r)
                if ~isempty(strfind(r(l2).name, 'white')) | ~isempty(strfind(r(l2).name, 'ventricle'))
                    thisv = spm_vol([NuisPath, r(l2).name]);
                    [thisY, thisXYZ] = spm_read_vols(thisv);
                    named(l).nuisroi(ncount).name = r(l2).name;
                    named(l).nuisroi(ncount).idx = find(thisY);
                    ncount = ncount + 1;
                else
                    disp('did not load any nuisance rois');
                end
            end
        end
    end
end

clear d r count fnames;     % Remove some unused variables

mask=1;

for l = 1:length(named) %cycle through the subjects
    tstart2=tic;   
    fprintf('\nSUBJECT %s\n',named(l).name);
    
    % Some subjects didn't use all 16 TRSE runs
    if isequal(project,'rest.nback') && isequal(task,'TRSE')
        if isequal(named(l).name,'106_A_high') || isequal(named(l).name,'107_A_low')
            blocklist = [5:16];
        end;
    end;
    
    if isequal(project,'task_graphtheory')
        if isequal(task,'rest')
            % Determine subject order from nback
            o=dir([EPIDIR,named(l).name,'/nback/afni_nback_*']);
            order=str2num(o(1).name(length(o(1).name)-1));
        else, % Either of the tasks
            % Determine subject order to get session list if nback or bimancoord
            o=dir([EPIDIR,named(l).name,EPIDIR_suff,'*']);
            for i=1:length(o),
                sessions{i}=o(i).name(length(o(i).name)-1);
            end;
            %%% Determine subject order based on session number
            if isequal(task,'nback')
                if isequal(sessions{1},'1'),
                    order=1;
                elseif isequal(sessions{1},'2'),
                    order=2;
                end;
            elseif isequal(task,'bimancoord')
                if isequal(sessions{1},'1'),
                    order=2;
                elseif isequal(sessions{1},'2'),
                    order=1;
                end;
            end;
        end;
    end;   
    
    % I needed to add a loop to go through multiple sessions
    for s = 1:length(sessions),
        disp(['    SUBJECT ' named(l).name ' SESSION ' sessions{s}]);
        sesscount=cellstr(num2str(s)); % Need this for despolarity file naming
        
        % Load mask to limit voxel analysis space
        disp(' ');
        if isequal(project,'Rest.Cerebellum')
            d = dir([EPIDIR, named(l).name, '/', EPIDIR_suff, 'r' named(l).name '-Mask.nii']);
        else
            if singleblock
                maskpref = [EPIDIR, named(l).name, EPIDIR_suff, num2str(blocklist) '/']; % This doesn't work for despolarity
            else
                if isequal(project,'despolarity')
                    maskpref = [EPIDIR, named(l).name, EPIDIR_suff, named(l).name, sessions{s} '/'];
                elseif isequal(project,'kki_adhd') && proc,
                    maskpref = [EPIDIR, tasks{1}, '/', named(l).fnum, '/']; % Can be from any task b/c same mask file is in all dirs
                elseif isequal(project,'task_graphtheory')
                    maskpref = [EPIDIR, named(l).name, maskdir_pref, named(l).name, maskdir_suff, normtype{s} '/'];
                elseif isequal(project,'rest.nback')
                    maskpref = [EPIDIR, named(l).name, EPIDIR_suff, '/'];
                else
                    maskpref = [EPIDIR, named(l).name, EPIDIR_suff '/'];
                end
            end
            if isequal(btwn_roi,'y')
                d = dir([maskpref, named(l).name '-Mask_resam.nii']);
            else
                if isequal(project,'despolarity') || isequal(project,'task_graphtheory')
                    d = dir([maskpref, named(l).name, '-Mask-ALLSESSRUNS.nii']); % Using mask that's a conjunction of ALL runs (across blocks and sessions); note this is seprate for A and B for task_graphtheory, as defined in maskpref above
                elseif isequal(project,'kki_adhd') && proc,
                    if task, % Rest and GNG
                        d = dir([maskpref,named(l).fnum,'_KKIproc_ns_mask_ALLRUNS.nii']);
                    else, % Rest only
                        d = dir([maskpref,named(l).fnum,'_RS_KKIproc_ns_mask.nii']);
                    end;
                elseif isequal(project,'kki_autism') && proc,
                    d = dir([maskpref,named(l).fnum,'.fsnwcfwa_ns_mask.nii']);
                elseif isequal(project,'rest.nback')
                    d = dir([maskpref, named(l).name, '-Mask-ALL',task,'.nii']); % Using mask that's a conjuntion of ALL runs per task
                else
                    d = dir([maskpref, named(l).name '-Mask.nii']);
                end
            end
        end
        meand = d;

        if ~isempty(d)
            fprintf('    Time elapsed: %.2f minutes for sub, %.2f minutes total.\n',toc(tstart2)/60,toc(tstart)/60);
            disp(['  Loading epi mask ', d(1).name]);
            % Add this back in if use a mean epi as a mask...but I'm using an already-binarized mask
            %disp(['    as a simple mask with ', num2str(100*maskthresh),'threshold...']);
            Vmask = spm_vol([maskpref,d(1).name]);
        else
            disp(['  Unable to load epi mask for subject ', named(l).name, '.']);
        end

        % Exclude voxels outside the brain
        [thisY, thisXYZ] = spm_read_vols(Vmask);
        if mask
            %mthresh = maskthresh * max(max(max(thisY)));
            %thisY(find(thisY <= mthresh)) = 0;
            named(l).maskvox = find(thisY);
            named(l).Ysize = size(thisY);
        else
            disp('Not masking against Mean EPI');
            named(l).maskvox = find(thisY);
            named(l).Ysize = size(thisY);
        end

        if isequal(btwn_roi,'y')
            d = dir([LesionDir, 'lesion_' named(l).name '.nii']);
            Vmask = spm_vol([LesionDir d.name]);
            [thisY, thisXYZ] = spm_read_vols(Vmask);
            named(l).lesmaskvox = find(thisY);
        end

        clear thisY thisXYZ;

        % Load filenames for EPIs, as well as the motion parameters
        fprintf('    Time elapsed: %.2f minutes for sub, %.2f minutes total.\n',toc(tstart2)/60,toc(tstart)/60);            
        disp(['  Loading epi files for subject ', named(l).name, '...']);
        
        epid=[];
        if isequal(project,'Atten.TMS') || isequal(project,'Atten.TMS2')
            if do_betacorr
                epid = dir([EPIDIR,named(l).name,EPIDIR_suff,'*' epi_name_1 named(l).name epi_name_2]); %this is for beta corrs
            else
                epid = dir([EPIDIR,named(l).name,EPIDIR_suff,'*' named(l).name '*' epi_name_1 epi_name_2]); %this is for whole blocks
            end
            blvals = [1:length(epid)];
            blocklist=[];
            for blval = 1:length(blvals)
                blocklist=[blocklist num2str(blvals(blval))];
            end
        elseif isequal(project,'Rest.Cerebellum')
            epid = dir([EPIDIR,named(l).name,EPIDIR_suff,'/r*' named(l).name '*' epi_name_1 epi_name_2]);
        elseif isequal(project,'Rest.LowTR')
            if singleblock
                for e = 1:length(blocklist)
                    epid = [epid; dir([EPIDIR,named(l).name,EPIDIR_suff,num2str(blocklist(e)),'/',named(l).name '*' epi_name_1 '*' blocklist(e) '*' epi_name_2])];
                end
            else
                for e = 1:length(blocklist)
                    epid = [epid; dir([EPIDIR,named(l).name,EPIDIR_suff,'/',named(l).name '*' epi_name_1 '*' blocklist(e) '*' epi_name_2])];
                end
            end
        else
            for e = 1:length(blocklist)
                % Need this loop b/c blocklist for these two subjects starts at 5, not 1, and need to account for the zero-padding correctly
                if isequal(project,'rest.nback') && (isequal(named(l).name,'106_A_high') || isequal(named(l).name,'107_A_low'))
                    double_digits=6;
                else
                    double_digits=10;
                end
                if e < double_digits
                    if isequal(project,'despolarity')
                        epid = [epid; dir([EPIDIR,named(l).name,EPIDIR_suff,named(l).name,sessions{s},'/',named(l).name '*' epi_name_1 num2str(blocklist(e)) epi_name_2])];
                    elseif isequal(project,'kki_adhd') && proc,
                        epid = [epid; dir([EPIDIR,tasks{e},'/',named(l).name,EPIDIR_suff{e},named(l).name,'_',tasktype{e},epi_name_2])];
                    elseif isequal(project,'task_graphtheory')
                        epid = [epid; dir([EPIDIR,named(l).name,EPIDIR_suff,named(l).name,sessions{s},normtype{s},'/',named(l).name '*' epi_name_1 num2str(blocklist(e)) blocksuff{e} epi_name_2])];
                    elseif isequal(project,'kki_autism') && proc,
                        epid = [epid; dir([EPIDIR,named(l).name,EPIDIR_suff,'/',named(l).name '*' epi_name_1 '*' epi_name_2])];
                    else
                        epid = [epid; dir([EPIDIR,named(l).name,EPIDIR_suff,'/',named(l).name '*' epi_name_1 '*0' num2str(blocklist(e)) '*' epi_name_2])];
                    end
                else 
                    epid = [epid; dir([EPIDIR,named(l).name,EPIDIR_suff,'/',named(l).name '*' epi_name_1 '*' num2str(blocklist(e)) '*' epi_name_2])];
                end
            end
        end
        %cd([EPIDIR,named(l).name,EPIDIR_suff]);
        %%%%%%%%%%%this for loop is new--do each block separate
        for epi=1:length(epid) %cycle through the epis for this subj
            
            if isequal(project,'rest.nback')
                disp(['    BLOCK ' num2str(blocklist(epi))]); % Because blocklist is double instead of string
            else
                disp(['    BLOCK ' blocklist(epi)]);
            end
            
            if singleblock
                EPIDIR_suff_new = [EPIDIR_suff,num2str(blocklist(epi)),'/'];
            else
                if isequal(project,'despolarity')
                    EPIDIR_suff_new = [EPIDIR_suff,named(l).name,sessions{s},'/'];
                elseif isequal(project,'task_graphtheory')
                    EPIDIR_suff_new = [EPIDIR_suff,named(l).name,sessions{s},normtype{s},'/'];
                else
                    EPIDIR_suff_new = [EPIDIR_suff,'/']; 
                end
            end
            if isequal(project,'kki_adhd')
                cd([EPIDIR,tasks{epi},'/',named(l).name,EPIDIR_suff_new{epi}]);
            else
                cd([EPIDIR,named(l).name,EPIDIR_suff_new]);
            end

            blocklength = length(spm_vol(epid(epi).name));
            if blocklength>1
                if isequal(project,'task_graphtheory')
                    actual_epi = num2str(epi); % So it's actual number instead of the name
                elseif isequal(project,'rest.nback')
                    actual_epi = num2str(blocklist(epi)); % So it's a string istead of a double
                else
                    actual_epi = blocklist(epi);
                end

                if isequal(project,'Atten.TMS') || isequal(project,'Atten.TMS2')
                    epitype = epid(epi).name(6:end-9); %this cuts off beta and subnum.nii
                end
                if isequal(project,'kki_adhd')
                    named(l).datafiles = [char(ones(length(meand(1)),1) * ...
                        [EPIDIR,tasks{epi},'/',named(l).name,EPIDIR_suff_new{epi}]),cat(1,epid(epi).name)];
                else
                    named(l).datafiles = [char(ones(length(meand(1)),1) * ...
                        [EPIDIR,named(l).name, EPIDIR_suff_new]),cat(1,epid(epi).name)];
                end

                if regress_nuis
                    if isequal(project,'Rest.LowTR')
                        d = dir([EPIDIR,named(l).name, EPIDIR_suff_new, motionfilename]);
                        named(l).motion = load([EPIDIR,named(l).name, EPIDIR_suff_new, d(1).name]);
                    elseif isequal(project,'task_graphtheory')
                        d = dir([EPIDIR,named(l).name, EPIDIR_suff_new, motionfilename,blocklist(epi),blocksuff{epi},'*',motionfilename_suff]);
                        named(l).motion = load([EPIDIR,named(l).name, EPIDIR_suff_new, d.name]);
                    else
                        d = dir([EPIDIR,named(l).name, EPIDIR_suff_new, motionfilename]);
                        if isequal(project,'Rest.TMSnew') || isequal(project,'MegaRest.TMS')
                           named(l).motion = load([EPIDIR,named(l).name, EPIDIR_suff_new, d(actual_epi).name]);
                        elseif isequal(project,'rest.nback')
                            named(l).motion = load([EPIDIR,named(l).name, EPIDIR_suff_new, d(epi).name]);
                        else
                           named(l).motion = load([EPIDIR,named(l).name, EPIDIR_suff_new, d(str2num(actual_epi)).name]);
                        end
                    end
                end
                % Load actual data
                Vdata = spm_vol(named(l).datafiles);
                [thisY, thisXYZ] = spm_read_vols(Vdata);
                for l3 = 1:size(thisY,4)
                    tmp = thisY(:,:,:,l3);
                    ALL(l).data(l3,:) = tmp(named(l).maskvox)';
                end

                clear tmp thisY thisXYZ;


                % Convert all ROI indices to the masked space indices
                fprintf('    Time elapsed: %.2f minutes for sub, %.2f minutes total.\n',toc(tstart2)/60,toc(tstart)/60);
                disp('converting rois');
                for l3 = 1:length(named(l).roi)
                    found = 0;
                    notfound = 0;
                    named(l).roi(l3).mask_space = [];
                    named(l).roi(l3).mask_missing = [];
                    for l4 = 1:length(named(l).roi(l3).idx)
                        tmpidx = find(named(l).maskvox == named(l).roi(l3).idx(l4));

                        %[I J] = ind2sub(named(l).maskvox == named(l).roi(l3).idx(l4));
                        if ~isempty(tmpidx)
                            found = found+1;
                            named(l).roi(l3).mask_space(found) = tmpidx;
                        else
                            notfound = notfound+1;
                            named(l).roi(l3).mask_missing(notfound) = named(l).roi(l3).idx(l4);
                        end
                    end
                    if isempty(named(l).roi(l3).mask_space)
                        disp(['     --> No voxels surpass threshold for ROI ', named(l).roi(l3).name]);
                    end
                    if ~isempty(named(l).roi(l3).mask_missing)
                        disp(['     --> Missing voxels [', int2str(named(l).roi(l3).mask_missing), ...
                            '] for ROI ', named(l).roi(l3).name]);
                        %if length(named(l).roi(l3).mask_missing)/length(named(l).roi(l3).mask_space) >= 0.75
                        if length(named(l).roi(l3).mask_missing)/length(named(l).roi(l3).idx) >= 0.75
                            named(l).roi(l3).mask_make = 0; %more than 75% of the ROI is outside the EPI
                            disp(['ROI ',named(l).roi(l3).name, ' missing more than 75% of voxels!!!!!!!!']);
                        else
                            named(l).roi(l3).mask_make = 1;
                        end
                    else
                        named(l).roi(l3).mask_make = 1;
                    end
                end

                if isequal(btwn_roi,'y') % Didn't change this loop for JRC projects
                    % Convert all ROI indices to the masked space indices
                    disp('converting within rois');
                    for l3 = 1:length(named(l).btwnroi)
                        found = 0;
                        notfound = 0;
                        found2 = 0; found3=0;
                        named(l).btwnroi(l3).mask_space = [];
                        named(l).btwnroi(l3).mask_missing = [];
                        named(l).btwnroi(l3).lesion_space = [];
                        named(l).btwnroi(l3).nonlesion_space = [];
                        named(l).btwnroi(l3).coords = [];
                        for l4 = 1:length(named(l).btwnroi(l3).idx)
                            tmpidx = find(named(l).maskvox == named(l).btwnroi(l3).idx(l4));
                            tmplesidx = find(named(l).lesmaskvox == named(l).btwnroi(l3).idx(l4));
                            named(l).btwnroi(l).coords(end+1,:) = Indices2coords(tmpidx,named(l).Ysize)';

                            if ~isempty(tmpidx)
                                found = found+1;
                                named(l).btwnroi(l3).mask_space(found) = tmpidx;
                                if ~isempty(tmplesidx)
                                    found2 = found2+1;
                                    named(l).btwnroi(l3).lesion_space(found2) = l4;%tmplesidx;
                                else
                                    found3 = found3+1;
                                    named(l).btwnroi(l3).nonlesion_space(found3) = l4;%tmpidx;
                                end
                            else
                                notfound = notfound+1;
                                named(l).btwnroi(l3).mask_missing(notfound) = named(l).btwnroi(l3).idx(l4);
                            end
                        end
                        if isempty(named(l).btwnroi(l3).mask_space)
                            disp(['     --> No voxels surpass threshold for ROI ', named(l).btwnroi(l3).name]);
                        end
                        if ~isempty(named(l).btwnroi(l3).mask_missing)
                            disp(['     --> Missing voxels [', int2str(named(l).btwnroi(l3).mask_missing), ...
                                '] for ROI ', named(l).btwnroi(l3).name]);
                            %if length(named(l).roi(l3).mask_missing)/length(named(l).roi(l3).mask_space) >= 0.75
                            if length(named(l).btwnroi(l3).mask_missing)/length(named(l).btwnroi(l3).idx) >= 0.75
                                named(l).btwnroi(l3).mask_make = 0; %more than 75% of the ROI is outside the EPI
                                disp(['ROI ',named(l).btwnroi(l3).name, ' missing more than 75% of voxels!!!!!!!!']);
                            else
                                named(l).btwnroi(l3).mask_make = 1;
                            end
                        else
                            named(l).btwnroi(l3).mask_make = 1;
                        end
                    end
                end

                if regress_nuis
                    for l3 = 1:length(named(l).nuisroi)
                        found = 0;
                        notfound = 0;
                        named(l).nuisroi(l3).mask_space = [];
                        named(l).nuisroi(l3).mask_missing = [];
                        for l4 = 1:length(named(l).nuisroi(l3).idx)
                            tmpidx = find(named(l).maskvox == named(l).nuisroi(l3).idx(l4));
                            if ~isempty(tmpidx)
                                found = found+1;
                                named(l).nuisroi(l3).mask_space(found) = tmpidx;
                            else
                                notfound = notfound+1;
                                named(l).nuisroi(l3).mask_missing(notfound) = named(l).nuisroi(l3).idx(l4);
                            end
                        end
                        if isempty(named(l).nuisroi(l3).mask_space)
                            disp(['     --> No voxels surpass threshold for ROI ', named(l).nuisroi(l3).name]);
                        end
                        if ~isempty(named(l).nuisroi(l3).mask_missing)
                            disp(['     --> Missing voxels [', int2str(named(l).nuisroi(l3).mask_missing), ...
                                '] for ROI ', named(l).nuisroi(l3).name]);
                        end
                    end
                end

                if do_betacorr || ((isequal(project,'kki_adhd') || isequal(project,'kki_autism')) && proc)
                    disp('No bandpass filter');
                else

                    % Apply bandpass filter
                    ALL(l).data=detrend(ALL(l).data);

                    disp('  Applying bandpass filter...');
                    for l3 = l:size(ALL(l).data,2)
                        ALL(l).data(:,l3) = fbandpass_new(ALL(l).data(:,l3), TR, f, 0);
                    end
                end

                if isequal(smooth,'y')
                    % Smooth filtered data
                    fprintf('    Time elapsed: %.2f minutes for sub, %.2f minutes total.\n',toc(tstart2)/60,toc(tstart)/60);
                    disp(['  Smoothing data with Gaussian of FWHM = ', int2str(p_smooth), '...']);
                    if ~exist([EPIDIR,named(l).name,EPIDIR_suff_new,'tmp_epis_' atlas '/'],'dir')
                        unix(['mkdir ', EPIDIR,named(l).name, EPIDIR_suff_new, 'tmp_epis_' atlas '/']);
                    end

                    tmp_epi_dir = [EPIDIR,named(l).name, EPIDIR_suff_new, 'tmp_epis_' atlas '/'];
                    results = Vdata(1);
                    for l3 = 1:size(ALL(l).data,1)
                        if l3 < 10
                            fnumber = ['000', int2str(l3)];
                        elseif l3 < 100
                            fnumber = ['00', int2str(l3)];
                        elseif l3 < 1000
                            fnumber = ['0', int2str(l3)];
                        else
                            fnumber = int2str(l3);
                        end
                        results.fname = [tmp_epi_dir, ['tmp_epis_' atlas], fnumber, '.nii'];%'.img'];
                        xxx = zeros(named(l).Ysize);
                        xxx(named(l).maskvox) = ALL(l).data(l3,:);

                        warning off all;
                        spm_write_vol(results, xxx);

                    end
                    for l3 = 1:size(ALL(l).data,1)
                        if l3 < 10
                            fnumber = ['000', int2str(l3)];
                        elseif l3 < 100
                            fnumber = ['00', int2str(l3)];
                        elseif l3 < 1000
                            fnumber = ['0', int2str(l3)];
                        else
                            fnumber = int2str(l3);
                        end

                        warning off all;
                        spm_smooth([tmp_epi_dir, ['tmp_epis_' atlas], fnumber, '.nii'],...
                            [tmp_epi_dir, 'stmp_epis', fnumber, '.nii'],p_smooth);
                    end
                    d = dir([tmp_epi_dir, 's*nii']);
                    simg = [char(ones(size(d,1),1) * tmp_epi_dir), cat(1,d.name)];
                    sVdata = spm_vol(simg);
                    [thisY, thisXYZ] = spm_read_vols(sVdata);
                    for l3 = 1:size(thisY,4)
                        tmp = thisY(:,:,:,l3);
                        ALL(l).data(l3,:) = tmp(named(l).maskvox)';
                    end
                    clear tmp thisY thisXYZ;
                    unix(['rm -R ', tmp_epi_dir]);
                else
                    disp('not smoothing');
                end

                % Get whole brain signal
                if wb_or_n == 'y'
                    fprintf('    Time elapsed: %.2f minutes for sub, %.2f minutes total.\n',toc(tstart2)/60,toc(tstart)/60);
                    disp('  Computing whole brain signal -- this MUST be better specified...');
                    named(l).wholebrain = mean(ALL(l).data,2);
                else
                    fprintf('    Time elapsed: %.2f minutes for sub, %.2f minutes total.\n',toc(tstart2)/60,toc(tstart)/60);
                    disp('  NOT including a whole brain signal as a nuisance regressor.');
                    named(l).wholebrain = [];
                end

                if regress_nuis
                    % Get the mean signal for the nuisance variables
                    for l3 = 1:length(named(l).nuisroi)
                        named(l).nuisroi(l3).data = mean(ALL(l).data(:,named(l).nuisroi(l3).mask_space),2);
                    end
                end

                tlen = size(ALL(l).data,1);

                % Diff_nuis approximates the temporal derivative
                if regress_nuis
                    % Begin regression to remove components related to nuisance variables
                    fprintf('    Time elapsed: %.2f minutes for sub, %.2f minutes total.\n',toc(tstart2)/60,toc(tstart)/60);
                    disp('  Beginning regression of nuisance variables...');
                    all_nuis = [named(l).motion, named(l).wholebrain];
                    for l3 = 1:length(named(l).nuisroi)
                        all_nuis = [all_nuis, named(l).nuisroi(l3).data];
                    end


                    diff_nuis = [zeros(1, size(all_nuis,2)); diff(all_nuis)];
                    all_nuis = [all_nuis, diff_nuis, ones(size(all_nuis,1),1)];

                    for l3 = 1:size(ALL(l).data,2)
                        b = regress(ALL(l).data(:,l3),all_nuis);
                        ALL(l).data(:,l3) = ALL(l).data(:,l3) - all_nuis*b;

                    end
                end
                % Get the mean signal for all ROIs
                for l3 = 1:length(named(l).roi)
                    if ~isempty(named(l).roi(l3).mask_space)
                        named(l).roi(l3).data = mean(ALL(l).data(:,named(l).roi(l3).mask_space),2);

                        %this is intended to save only the data within the ROI
                        %(before taking the mean)
                        named(l).roi(l3).data2 = ALL(l).data(:,named(l).roi(l3).mask_space);

                        if cf_or_n == 'c'
                            named(l).roi(l3).correction = sum(xcorr( (named(l).roi(l3).data - ...
                                mean(named(l).roi(l3).data)), 'coeff').^2);
                        else
                            named(l).roi(l3).correction = mias_cf;
                        end
                    else
                        named(l).roi(l3).data = zeros(size(mean(ALL(l).data(:,named(l).roi(l3).mask_space),2)));
                    end
                end

                if ~do_betacorr & isequal(atlas,'dosenbach') % This loop was for all atlases, but is Dosenbach-specific so I added that clause here (JRC)
                    % Find the average signal across Dosenbach nodes (except for the ROI of interest)
                    count = 0;
                    for l3 = 1:length(named(l).roi)

                        if l3 < 8
                            tmpdata = [];
                            %CO network
                            for l4 = 1:7
                                if l3 ~=l4
                                    tmpdata = [tmpdata named(l).roi(l4).data];
                                end
                            end
                            COdata(:,l3) = mean(tmpdata,2);
                        else
                            tmpdata = [];
                            %FP network
                            for l4 = 8:18

                                if l3 ~=l4
                                    tmpdata = [tmpdata named(l).roi(l4).data];
                                end
                            end
                            count = count+1;
                            FPdata(:,count) = mean(tmpdata,2);
                        end
                    end
                end

                %%%%%%START RUNNING CORRELATIONS%%%%%%
                blocklength_orig = blocklength;
                if isequal(project,'Rest.TMSnew') || isequal(project,'MegaRest.TMS') % actual_epi is a double for these studies
                    startepi = (actual_epi-1)*(reps-1);
                else
                    startepi = (str2num(actual_epi)-1)*(reps-1);
                end
                for l5=1:reps
                    % Note that don't need these special loops if I define everything generally as I did below
                    %if reps==2
                    %    if isequal(project,'Rest.TMSnew') || isequal(project,'MegaRest.TMS')
                    %        blocknum=actual_epi*2 - mod(l5,2);
                    %    else
                    %        blocknum=str2num(actual_epi)*2 - mod(l5,2);
                    %    end
                    %    blocklength=floor(blocklength_orig/2);
                    %elseif reps == 4 || reps ==5
                    %    if isequal(project,'Rest.TMSnew') || isequal(project,'MegaRest.TMS')
                    %        blocknum=actual_epi + l5-1 + startepi;
                    %    else
                    %        blocknum=str2num(actual_epi) + l5-1 + startepi;
                    %    end
                    %    blocklength=floor(blocklength_orig/reps);
                    %else % Note that this assumes reps==1, when it can be any other value (other than 2,4, or 5) and go into this loop as well
                    %    if isequal(project,'Rest.TMSnew') || isequal(project,'MegaRest.TMS')
                    %        blocknum=actual_epi;
                    %    else
                    %        blocknum=str2num(actual_epi) + l5-1 + startepi;
                    %    end
                    %end
                    if isequal(project,'Rest.TMSnew') || isequal(project,'MegaRest.TMS') || isequal(project,'MegaRest.TMS_ASLconn')
                        blocknum=actual_epi;
                    else
                        blocknum=str2num(actual_epi) + l5-1 + startepi;
                    end
                    blocklength=floor(blocklength_orig/reps);

                    if l5==1
                        start=1;
                        stop=blocklength;
                    else
                        start = stop + 1;
                        stop = start + blocklength-1;
                    end

                    if ~isequal(project,'Atten.TMS') && ~isequal(project,'Atten.TMS2')
                        if isequal(project,'Rest.thetaTMS')
                            if isequal(prepost,'post')
                                blocknum = blocknum+8;
                            end
                        end
                    end
                    if blocknum < 10
                        strblock=['0' num2str(blocknum)];
                    else
                        strblock=num2str(blocknum);
                    end
                    epitype = ['Block' strblock];

                    if l5 >= 4 && isequal(project,'Rest.TMS') && isequal(named(l).name,'101') && isequal(actual_epi,'4')
                        disp('skipping these');
                    else
                        if isequal(vox2vox,'n') || isequal(vox2vox,'b')
                            % Compute the correlations for each of the ROIs with all the other ROIs
                            fprintf('    Time elapsed: %.2f minutes for sub, %.2f minutes total.\n',toc(tstart2)/60,toc(tstart)/60);
                            disp('  Computing correlations across all ROIs...');

                            data=[];

                            if isequal(win_roi,'y') % This loop has not been adapted for JRC projects
                                %within ROI correlations
                                for l3=1:length(named(l).roi)
                                    data=[];
                                    temproidir = [StudyPath, StudyPathnm, 'corr_'  named(l).roi(l3).name(6:end-8) '_mask/'];
                                    if ~exist(temproidir,'dir')
                                        unix(['mkdir ', temproidir]);
                                    end
                                    tempfilename = [temproidir, named(l).name, '_' epitype];
                                    if ~exist([tempfilename '.mat'],'file')
                                        disp(['Working on: ', tempfilename ]);
                                        for l4=1:size(named(l).roi(l3).mask_space,2)
                                            if named(l).roi(l3).mask_make
                                                data = [data ALL(l).data(start:stop,named(l).roi(l3).mask_space(l4))];
                                            else
                                                data = [data zeros(size(ALL(l).data(start:stop,named(l).roi(l3).mask_space(l4))))];
                                            end
                                        end
                                        [r] = corr(data);
                                        fdata = r;
                                        save(tempfilename, 'fdata');
                                        clear r fdata;
                                    else
                                        disp('file exists');
                                    end
                                end
                            elseif isequal(btwn_roi,'y') % This loop has not been adapted for JRC projects
                                for rval = 1:length(named(l).btwnroi.idx)
                                    if rval <10
                                        strrval = ['00' num2str(rval)];
                                    elseif rval <100
                                        strrval = ['0' num2str(rval)];
                                    else
                                        strrval = num2str(rval);
                                    end
                                    tempfilename = [StudyPath, StudyPathnm, 'vox2roi/', named(l).fnum, '_' epitype, '/' atlas btwn_roi_val '/vox',strrval];
                                    %data=[ALL(l).data(start:stop,named(l).btwnroi.mask_space(rval))];
                                    data = [];
                                    if ~exist(tempfilename)
                                        for l3=1:length(named(l).roi)
                                            if named(l).roi(l3).mask_make
                                                if l3 == str2num(btwn_roi_val)
                                                    data=[data ALL(l).data(start:stop,named(l).btwnroi.mask_space(rval))];
                                                else
                                                    data=[data named(l).roi(l3).data(start:stop)];
                                                end
                                            else
                                                data=[data zeros(size(named(l).roi(l3).data(start:stop)))];
                                            end
                                        end

                                        [r] = corr(data);
                                        fdata = r;


                                        save(tempfilename, 'fdata');
                                        clear r fdata;
                                    else
                                        disp(['exists: ', tempfilename]);
                                    end
                                end
                            else
                                if do_betacorr 
                                    tempfilename = [AnalPath, StudyPathnm, 'roi2roi/', named(l).fnum, '_' beta_cond];
                                else
                                    if isequal(project,'task_graphtheory')
                                        % Switch rest and task session numbers that will be in the file names so
                                        % they're the same as if order was 1; don't change them if order == 1. 
                                        % This is so numbers are consistent when I want to relate rest and task 
                                        % runs within a session, so all nback days are labeled 1 and all biman 
                                        % days are labeled 2, even if it was really 2/1. Don't change real session
                                        % session number though so rest of script will run correctly!!
                                        if order == 2, 
                                            if isequal(sessions{s},'1'),
                                                sessname='2';
                                            elseif isequal(sessions{s},'2'),
                                                sessname='1';
                                            end;
                                        elseif order == 1,
                                            if isequal(sessions{s},'1'),
                                                sessname='1';
                                            elseif isequal(sessions{s},'2'),
                                                sessname='2';
                                            end;
                                        end;
                                    end;
                                    if isequal(smooth,'y') || proc, % Already processed so smoothed before this script
                                        if isequal(project,'despolarity') || isequal(project,'kki_adhd') || isequal(project,'task_graphtheory') || isequal(project,'kki_autism') || isequal(project,'rest.nback')
                                            % Make directories that may not exist
                                            if ~exist([AnalPath])
                                                system(['mkdir ' AnalPath]);
                                            end
                                            if ~exist([AnalPath, StudyPathnm])
                                                system(['mkdir ' AnalPath,StudyPathnm]);
                                            end
                                            if ~exist([AnalPath, StudyPathnm, 'roi2roi/'])
                                                system(['mkdir ' AnalPath,StudyPathnm, 'roi2roi/']);
                                            end
                                            if isequal(project,'despolarity') % Use sesscount instead of sessions{s} to account for some subs who I'm not using consecutive sessions for
                                                tempfilename = [AnalPath, StudyPathnm, 'roi2roi/', named(l).fnum, '_Sess' sesscount{1} '_' epitype];
                                                temptxtfilename = [AnalPath, StudyPathnm, 'roi2roi/TXTfiles/', named(l).fnum, '_Sess' sesscount{1} '_' epitype,'.txt'];
                                                tempTSfilename = [AnalPath, StudyPathnm, 'roi2roi/timeseries/', named(l).fnum, '_Sess' sesscount{1} '_' epitype '_timeseries'];
                                                tempTStxtfilename = [AnalPath, StudyPathnm, 'roi2roi/timeseries/TXTfiles/', named(l).fnum, '_Sess' sesscount{1} '_' epitype,'_timeseries.txt'];
                                            elseif isequal(project,'task_graphtheory') % Use sessname instead of sessions{s} so rest can be associated with appropriate task day
                                                tempfilename = [AnalPath, StudyPathnm, 'roi2roi/', task, '_', normtype{s}, '/', named(l).fnum, '_Sess' sessname '_' epitype];
                                                temptxtfilename = [AnalPath, StudyPathnm, 'roi2roi/', task, '_', normtype{s}, '/TXTfiles/', named(l).fnum, '_Sess' sessname '_' epitype,'.txt'];
                                                tempTSfilename = [AnalPath, StudyPathnm, 'roi2roi/', task, '_', normtype{s}, '/timeseries/', named(l).fnum, '_Sess' sessname '_' epitype '_timeseries'];
                                                tempTStxtfilename = [AnalPath, StudyPathnm, 'roi2roi/', task, '_', normtype{s}, '/timeseries/TXTfiles/', named(l).fnum, '_Sess' sessname '_' epitype '_timeseries.txt'];
                                            elseif isequal(project,'rest.nback')
                                                tempfilename = [AnalPath, StudyPathnm, 'roi2roi/', task, '/', named(l).fnum, '_Sess' sessions{s} '_' epitype];
                                                temptxtfilename = [AnalPath, StudyPathnm, 'roi2roi/', task, '/TXTfiles/', named(l).fnum, '_Sess' sessions{s} '_' epitype,'.txt'];
                                                tempTSfilename = [AnalPath, StudyPathnm, 'roi2roi/', task, '/timeseries/', named(l).fnum, '_Sess' sessions{s} '_' epitype '_timeseries'];
                                                tempTStxtfilename = [AnalPath, StudyPathnm, 'roi2roi/', task, '/timeseries/TXTfiles/', named(l).fnum, '_Sess' sessions{s} '_' epitype,'_timeseries.txt'];
                                            else
                                                tempfilename = [AnalPath, StudyPathnm, 'roi2roi/', named(l).fnum, '_Sess' sessions{s} '_' epitype];
                                                temptxtfilename = [AnalPath, StudyPathnm, 'roi2roi/TXTfiles/', named(l).fnum, '_Sess' sessions{s} '_' epitype,'.txt'];
                                                tempTSfilename = [AnalPath, StudyPathnm, 'roi2roi/timeseries/', named(l).fnum, '_Sess' sessions{s} '_' epitype '_timeseries'];
                                                tempTStxtfilename = [AnalPath, StudyPathnm, 'roi2roi/timeseries/TXTfiles/', named(l).fnum, '_Sess' sessions{s} '_' epitype,'_timeseries.txt'];
                                            end
                                        else
                                            tempfilename = [StudyPath, StudyPathnm, 'roi2roi/', named(l).fnum, '_' epitype];
                                            temptxtfilename = [StudyPath, StudyPathnm, 'roi2roi/TXTfiles/', named(l).fnum, '_' epitype,'.txt'];
                                        end
                                    else
                                        if isequal(project,'despolarity'),
                                            tempfilename = [AnalPath, StudyPathnm, 'roi2roi_unsmooth/', named(l).fnum, '_Sess' sesscount{1} '_' epitype];
                                            temptxtfilename = [AnalPath, StudyPathnm, 'roi2roi_unsmooth/TXTfiles/', named(l).fnum, '_Sess' sesscount{1} '_' epitype,'.txt'];
                                        elseif isequal(project,'kki_adhd') || isequal(project,'kki_autism'),
                                            tempfilename = [AnalPath, StudyPathnm, 'roi2roi_unsmooth/', named(l).fnum, '_Sess' sessions{s} '_' epitype];
                                            temptxtfilename = [AnalPath, StudyPathnm, 'roi2roi_unsmooth/TXTfiles/', named(l).fnum, '_Sess' sessions{s} '_' epitype,'.txt'];
                                        elseif isequal(project,'task_graphtheory') % Use sessname instead of sessions{s} so rest can be associated with appropriate task day
                                            tempfilename = [AnalPath, StudyPathnm, 'roi2roi_unsmooth/', task, '_', normtype{s}, '/', named(l).fnum, '_Sess' sessname '_' epitype];
                                            temptxtfilename = [AnalPath, StudyPathnm, 'roi2roi_unsmooth/', task, '_', normtype{s}, '/TXTfiles/', named(l).fnum, '_Sess' sessname '_' epitype,'.txt'];
                                        elseif isequal(project,'rest.nback')
                                            tempfilename = [AnalPath, StudyPathnm, 'roi2roi_unsmooth/', task, '/', named(l).fnum, '_Sess' sessions{s} '_' epitype];
                                            temptxtfilename = [AnalPath, StudyPathnm, 'roi2roi_unsmooth/', task, '/TXTfiles/', named(l).fnum, '_Sess' sessions{s} '_' epitype,'.txt'];
                                        else
                                            tempfilename = [StudyPath, StudyPathnm, 'roi2roi_unsmooth/', named(l).fnum, '_' epitype];
                                        end
                                    end
                                end
                                if ~exist(tempfilename) || ~exist(tempTSfilename)
                                    for l3=1:length(named(l).roi)
                                        if named(l).roi(l3).mask_make
                                            data=[data named(l).roi(l3).data(start:stop)];
                                        else
                                            data=[data zeros(size(named(l).roi(l3).data(start:stop)))];
                                        end
                                    end

                                    [r] = corr(data);
                                    fdata = r;
                                    
                                    %%%% FOR NOW ONLY DOING r BECAUSE THAT'S WHAT IS NEEDED FOR GRAPH THEORY...CAN INCORPORATE THIS LATER
                                    % What Andy used..equivalent to atanh(r)*sqrt(N-3)/CF-not sure why he
                                    % divided by CF; I'm going to NOT divide by CF (Emi doesn't)
                                    %zscore = 0.5 * log ((1+r)./(1-r));
                                    %zrho = (1/sqrt((tlen-3)/named(l).roi(l3).correction));
                                    %Z = zscore./zrho;
                                    Z = atanh(r)*sqrt(tlen-3);
                                    %%%%
                                    
                                    if isequal(project,'despolarity') || isequal(project,'kki_adhd') || isequal(project,'kki_autism') 
                                        temptxtfiledir = [AnalPath, StudyPathnm, 'roi2roi/TXTfiles/'];
                                        tempTSfiledir = [AnalPath, StudyPathnm, 'roi2roi/timeseries/'];
                                        tempTStxtfiledir = [AnalPath, StudyPathnm, 'roi2roi/timeseries/TXTfiles/'];
                                    elseif isequal(project,'task_graphtheory')
                                        tempfiledir = [AnalPath, StudyPathnm, 'roi2roi/' task '_' normtype{s} '/'];
                                        temptxtfiledir = [AnalPath, StudyPathnm, 'roi2roi/' task '_' normtype{s} '/TXTfiles/'];
                                        tempTSfiledir = [AnalPath, StudyPathnm, 'roi2roi/' task '_' normtype{s} '/timeseries/'];
                                        tempTStxtfiledir = [AnalPath, StudyPathnm, 'roi2roi/' task '_' normtype{s} '/timeseries/TXTfiles/'];
                                    elseif isequal(project,'rest.nback')
                                        tempfiledir = [AnalPath, StudyPathnm, 'roi2roi/' task '/'];
                                        temptxtfiledir = [AnalPath, StudyPathnm, 'roi2roi/' task '/TXTfiles/'];
                                        tempTSfiledir = [AnalPath, StudyPathnm, 'roi2roi/' task '/timeseries/'];
                                        tempTStxtfiledir = [AnalPath, StudyPathnm, 'roi2roi/' task '/timeseries/TXTfiles/'];
                                    else
                                        temptxtfiledir = [StudyPath, StudyPathnm, 'roi2roi/TXTfiles/'];
                                        tempTSfiledir = [AnalPath, StudyPathnm, 'roi2roi/timeseries/'];
                                        tempTStxtfiledir = [AnalPath, StudyPathnm, 'roi2roi/timeseries/TXTfiles/'];
                                    end
                                    if isequal(project,'task_graphtheory') || isequal(project,'rest.nback')
                                        if ~exist(tempfiledir)
                                            system(['mkdir ' tempfiledir]);
                                        end
                                    end
                                    if ~exist(temptxtfiledir)
                                        system(['mkdir ' temptxtfiledir]);
                                    end
                                    if ~exist(tempTSfiledir)
                                        system(['mkdir ' tempTSfiledir]);
                                    end
                                    if ~exist(tempTStxtfiledir)
                                        system(['mkdir ' tempTStxtfiledir]);
                                    end
                                    % Save timeseries
                                    save(tempTSfilename, 'data');
                                    save(tempTStxtfilename, 'data','-ascii');
                                    % Save correlation matrices
                                    save(tempfilename, 'fdata');
                                    save(temptxtfilename,'fdata','-ascii');
                                    clear r fdata;
                                else
                                    disp(['exists: ', tempfilename]);
                                end
                                if isequal(partcorr,'y') % This loop has not been adapted for JRC projects
                                    if 0
                                        for l4=1:length(named(l).roi)
                                            partdata = named(l).roi(l4).data(start:stop);
                                            data=[];
                                            for l3=1:length(named(l).roi)
                                                if named(l).roi(l3).mask_make
                                                    data=[data named(l).roi(l3).data(start:stop)];
                                                else
                                                    data=[data zeros(size(named(l).roi(l3).data(start:stop)))];
                                                end
                                            end

                                            [r] = partialcorr(data,partdata);
                                            fdata = r;
                                            tempfilename = [StudyPath, StudyPathnm, 'roi2roi/part_', named(l).roi(l4).name(1:end-8),'_',named(l).fnum, '_' epitype];
                                            save(tempfilename, 'fdata');
                                            clear r fdata;
                                        end
                                    end

                                    %regress out all rois
                                    outmat = zeros(length(named(l).roi),length(named(l).roi));

                                    % make an array with all ROI time series
                                    alldatavals = [];
                                    for l3=1:length(named(l).roi)
                                        alldatavals = [alldatavals named(l).roi(l3).data(start:stop)];
                                    end
                                    allroilist = [1:length(named(l).roi)];
                                    for l3=1:length(named(l).roi)
                                        for l4=1:length(named(l).roi)
                                            npair1= allroilist(allroilist~=l3);
                                            npair2= allroilist(allroilist~=l4);
                                            nonpairind = intersect(npair1,npair2);

                                            partdata = alldatavals(:,nonpairind);

                                            if named(l).roi(l3).mask_make && named(l).roi(l4).mask_make
                                                data=[named(l).roi(l3).data(start:stop) named(l).roi(l4).data(start:stop)];
                                            else
                                                data=[zeros(size(named(l).roi(l3).data(start:stop))) zeros(size(named(l).roi(l3).data(start:stop)))];
                                            end

                                            outval = partialcorr(data,partdata);
                                            outmat(l3,l4) = outval(1,2);
                                        end
                                    end

                                    fdata = outmat;
                                    tempfilename = [StudyPath, StudyPathnm, 'roi2roi/partall_',named(l).name, '_' epitype];
                                    save(tempfilename, 'fdata');
                                    clear outmat fdata;
                                end
                            end
                        end

                        if isequal(vox2vox,'y') || isequal(vox2vox,'b') % This loop has not been adapted for JRC projects

                            ALL_old = ALL; %save this so do not overwrite
                            named_old = named; %save this so do not overwrite

                            % Compute the correlations for each of the ROIs across all
                            % voxels within the ROIs
                            disp('  Computing correlations across all voxels...');

                            if isequal(mask_epi_data,'y')

                                %this allows you to only include a portion of the brain
                                for mask_dfile = 1:size(mask_epi_data_file,1)

                                    for hemi = 1:size(hemi_list,1)
                                        disp(['Mask: ', mask_epi_data_file(mask_dfile,:), ', Hemi: ', hemi_list(hemi,:)]);
                                        tempMaskfilename = [MaskPath, hemi_list(hemi,:), '_', mask_epi_data_file(mask_dfile,:), '_', named(l).name, '.nii'];
                                        Vmask = spm_vol(tempMaskfilename);

                                        [thisY, thisXYZ] = spm_read_vols(Vmask);
                                        maskvox = find(thisY);

                                        ALL(l).data = ALL_old(l).data(:,maskvox);
                                        named(l).Ysize = size(thisY);
                                        named(l).maskvox = maskvox;

                                        %iterate through ROIs
                                        for l3 = 1:length(named(l).roi)
                                            if ~isempty(named(l).roi(l3).mask_space)
                                                %disp(['  Correction factor for ROI ', named(l).roi(l3).name, ' = ', ...
                                                %    num2str(named(l).roi(l3).correction)]);

                                                for l4 = 1:size(ALL(l).data,2)
                                                    r(l4) = corr(named(l).roi(l3).data(start:stop),ALL(l).data(start:stop,l4));
                                                end

                                                % What Andy used..equivalent to atanh(r)*sqrt(N-3)/CF-not sure why he
                                                % divided by CF; I'm going to NOT divide by CF (Emi doesn't)
                                                %zscore = 0.5 * log ((1+r)./(1-r));
                                                %zrho = (1/sqrt((tlen-3)/named(l).roi(l3).correction));
                                                %Z = zscore./zrho;
                                                Z = atanh(r)*sqrt(tlen-3);
                                                
                                                x = zeros(named(l).Ysize);
                                                xorig = x;
                                                xorig(named(l).maskvox) = r'; %CG
                                                x(named(l).maskvox) = Z'; %CG


                                                results = Vdata(1);
                                                results2=Vdata(1);

                                                if ~exist([StudyPath, StudyPathnm],'dir');
                                                    unix(['mkdir ', StudyPath, StudyPathnm]);
                                                end

                                                results.fname = [StudyPath, StudyPathnm, 'roi2vox/',hemi_list(hemi,:),'_', mask_epi_data_file(mask_dfile,:) '_' epitype, ...
                                                    '_', named(l).roi(l3).name];
                                                if isequal(output,'z')
                                                    dummy = spm_write_vol(results, x);
                                                elseif isequal(output,'r')
                                                    dummy = spm_write_vol(results,xorig);
                                                elseif isequal(output,'b')
                                                    dummy = spm_write_vol(results, x);
                                                    dummy2 = spm_write_vol(results2,xorig);
                                                end

                                                clear b r x tmpY tmpXYZ newv;
                                            end
                                        end
                                    end
                                end
                            
                            else
                                %computing on the whole brain
                                disp('roi to voxel!');
                                for l3 = 1:length(named(l).roi)
                                    if ~isempty(named(l).roi(l3).mask_space)
                                        if do_betacorr
                                            fileoutname = [StudyPath, StudyPathnm, 'roi2vox/WB_', beta_cond, '_',named(l).roi(l3).name];
                                        else
                                            fileoutname = [StudyPath, StudyPathnm, 'roi2vox/WB_', ...
                                                epitype, '_', named(l).roi(l3).name];
                                        end
                                        if exist(fileoutname,'file') || exist([fileoutname '.gz'],'file')
                                            disp(['file exists: ', fileoutname]);
                                        else
                                            disp(['making: ',fileoutname]);
                                            for l4 = 1:size(ALL(l).data,2)
                                                r(l4) = corr(named(l).roi(l3).data(start:stop),ALL(l).data(start:stop,l4));
                                            end

                                            % What Andy used..equivalent to atanh(r)*sqrt(N-3)/CF-not sure why he
                                            % divided by CF; I'm going to NOT divide by CF (Emi doesn't)
                                            %zscore = 0.5 * log ((1+r)./(1-r));
                                            %zrho = (1/sqrt((tlen-3)/named(l).roi(l3).correction));
                                            %Z = zscore./zrho;
                                            Z = atanh(r)*sqrt(tlen-3);

                                            x = zeros(named(l).Ysize);
                                            xorig = x;
                                            xorig(named(l).maskvox) = r'; %CG
                                            x(named(l).maskvox) = Z'; %CG


                                            results = Vdata(1);
                                            results2=Vdata(1);

                                            if ~exist([StudyPath, StudyPathnm],'dir');
                                                unix(['mkdir ', StudyPath, StudyPathnm]);
                                            end

                                            results.fname = fileoutname;
                                            if isequal(output,'z')
                                                dummy = spm_write_vol(results, x);
                                            elseif isequal(output,'r')
                                                dummy = spm_write_vol(results,xorig);
                                            elseif isequal(output,'b')
                                                dummy = spm_write_vol(results, x);
                                                dummy2 = spm_write_vol(results2,xorig);
                                            end

                                            clear b r x tmpY tmpXYZ newv;
                                        end
                                    end
                                end
                                if dodosen
                                    disp('computing dosenbach average maps');
                                    for l3 = 1:length(named(l).roi)
                                        if ~isempty(named(l).roi(l3).mask_space)
                                            fileoutname = [StudyPath, StudyPathnm, 'roi2vox/Dosen_WB_', ...
                                                epitype, '_', named(l).roi(l3).name];
                                            if exist(fileoutname,'file') || exist([fileoutname '.gz'],'file')
                                                disp(['file exists: ', fileoutname]);
                                            else
                                                disp(['making: ',fileoutname]);
                                                for l4 = 1:size(ALL(l).data,2)
                                                    if l3 < 8
                                                        r(l4) = corr(COdata(:,l3),ALL(l).data(start:stop,l4));
                                                    else
                                                        r(l4) = corr(FPdata(:,l3-7),ALL(l).data(start:stop,l4));
                                                    end
                                                    %r(l4) = corr(named(l).roi(l3).data(start:stop),ALL(l).data(start:stop,l4));
                                                end

                                                % What Andy used..equivalent to atanh(r)*sqrt(N-3)/CF-not sure why he
                                                % divided by CF; I'm going to NOT divide by CF (Emi doesn't)
                                                %zscore = 0.5 * log ((1+r)./(1-r));
                                                %zrho = (1/sqrt((tlen-3)/named(l).roi(l3).correction));
                                                %Z = zscore./zrho;
                                                Z = atanh(r)*sqrt(tlen-3);
                        
                                                x = zeros(named(l).Ysize);
                                                xorig = x;
                                                xorig(named(l).maskvox) = r'; %CG
                                                x(named(l).maskvox) = Z'; %CG


                                                results = Vdata(1);
                                                results2=Vdata(1);

                                                if ~exist([StudyPath, StudyPathnm],'dir');
                                                    unix(['mkdir ', StudyPath, StudyPathnm]);
                                                end

                                                results.fname = fileoutname;
                                                if isequal(output,'z')
                                                    dummy = spm_write_vol(results, x);
                                                elseif isequal(output,'r')
                                                    dummy = spm_write_vol(results,xorig);
                                                elseif isequal(output,'b')
                                                    dummy = spm_write_vol(results, x);
                                                    dummy2 = spm_write_vol(results2,xorig);
                                                end

                                                clear b r x tmpY tmpXYZ newv;
                                            end
                                        end
                                    end
                                end
                            end

                            ALL = ALL_old; %save this so do not overwrite
                            named = named_old; %save this so do not overwrite
                        % These lines are unecessary, because they're printed if I don't want to do vox2vox, but theroi2roi loop was already completed
                        %else
                        %    disp('what did you mean to do?');
                        end

                    end
                end
                clear ALL;
                if isequal(btwn_roi,'y') % This loop has not been adapted for JRC projects
                    save([StudyPath, StudyPathnm, 'vox2roi/', named(l).fnum, '_' epitype, '/' atlas btwn_roi_val '/named'],'named');
                end
                
                fprintf('    Time elapsed: %.2f minutes for sub, %.2f minutes total.\n',toc(tstart2)/60,toc(tstart)/60);
                    
            else
                disp('not enough trials here'); % this goes with the if blocklength>1 loop...so saying not enough trials if only 1 acquisition
            end
        end % end epi loop
    end % end session loop        
end % end sub loop

cd([ScriptPath]);
