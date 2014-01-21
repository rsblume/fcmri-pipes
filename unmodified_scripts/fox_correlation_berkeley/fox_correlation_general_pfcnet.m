function [] = fox_correlation_general_new(atlas,project_type,partcorr)

% The method is taken from Fox et.al., PNAS, 7/5/2005
% The steps are as follows:
%
% 1) Load EPI data.
% 2) Band-pass filter the time series (0.009 < f < 0.08)
% 3) Load nuisance parameters
%    a) rigid-body motion parameters
%    b) whole-brain mean signal, if desired
%    c) ventricular signal
%    d) white matter signal
% 4) Perform linear regression to remove the above -- i.e. y = Bx, where y
%       is the time series for a given voxel, x is the matrix of
%       regressors, and E = y - Bx is the desired output (I believe).
% 5) Do correlation maps across the whole brain, using a seed region taken
%       from previous contrasts (e.g. R frontal ROI).
% 6) Determine the degrees of freedom = #time points/correction factor.
%       Correction factor from Jenkins & Watts, Spectral Analysis & Its
%       Applications, 1968 (Holden-Day, San Francisco).
% 7) Calculate z-scores, with std.dev = 1/sqrt(df - 3)
%
% Inputs: 1) atlas: ROI folder
%         2) project_type: 'A', 'B', etc. depending on the study (from
%         params.txt)
%
% Outputs: None to calling function

% FOR SOME REASON, DEFAULTS ARE NOT LOADED BY SPM_DEFAULT CALL.
% *************************************************
%%%%%%%%%%%%%%%%%%%%%%
spm_defaults;
p_smooth = 6;
f = [0.009 0.08];
maskthresh = .1;
mias_cf = 2.34;     % Was Mia's value of 6, but concerned it's changed.  Will
%   use Fox's value for now, but needs to be recomputed.

disp(' ');

%%%%%%%%%%%%%%%%%%%%

allsubj = 'n';
win_roi = 'n';
vox2vox = 'n';

output = 'r'; %or 'z' or 'b'



if isequal(partcorr,1)
    ROIgroups = {''};
    hemis = {''};
elseif isequal(partcorr,2)
    ROIgroups = {''};
    hemis = {'_L'};
elseif isequal(partcorr,3)
    ROIgroups = {'DM','FDL','FVL'};%{'DL','DM','DRL','FDL','FVL','VL'};
    hemis = {'_L'};
    
elseif isequal(partcorr,4)
    ROIgroups = {'DM','FDL','FVL'};%{'DL','DM','DRL','FDL','FVL','VL'};
    hemis = {'_L'};
    
elseif isequal(partcorr,5)
    ROIgroups = {'DM','FDL','FVL'};%{'DL','DM','DRL','FDL','FVL','VL'};
    hemis = {'_L'};
    
else
    disp('no partial correlation');
    ROIgroups = {''};
    hemis = {''};
end

%use these for partial level correlations
%FDL_Lind = [47:2:57];
%FVL_Lind = [59:2:69];
%DM_Lind = [13:2:23];

if isequal(project_type,'sham')
    StudyPath = ['/home/despo/rstate/data/Rest.ShamTMS/Data/'];
    NuisPath = [StudyPath 'Masks/nuisance/masked/'];
    EPIDIR = StudyPath;
    EPIDIR_suff = '/Total/NIfTI/';
    motionfilename = '*00*-1D.txt';
    epi_name_1 = 'EPI-00*';
    epi_name_2 = '*-CoReg.nii';
    regress_nuis = 1;
    
    TR = 1.37;
    subjects = [114,116:119,201,203:220];
    reps = 2;
else
    StudyPath = ['/home/despo/rstate/data/Rest.PFCnet/Data/'];
    NuisPath = [StudyPath 'Masks/nuisance/'];
    EPIDIR = StudyPath;
    EPIDIR_suff = '/Analysis/Rest/afni_rest/';
    motionfilename = 'dfile.r*.1D';
    epi_name_1 = '*EPI-00';
    epi_name_2 = '*.nii';
    regress_nuis = 1;
    
    if isequal(project_type,'A')
        TR = 1;
        subjects = [401:414];
        reps = 1;
    elseif isequal(project_type,'B')
        TR = 2;
        reps = 1;
        subjects = [415:420,422:438];
    elseif isequal(project_type,'C')
        TR = 2;
        subjects = [439:445,448:459];
        reps = 1;
    elseif isequal(project_type,'D')
        TR = 1.37;
        subjects = [460:465,493:499];
        reps = 1;
    elseif isequal(project_type,'E')
        TR = 2;
        subjects = [466:484];
        reps = 1;
    elseif isequal(project_type,'F')
        TR = 1.37;
        subjects = [485];
        reps = 2;
    elseif isequal(project_type,'G')
        TR = 1.37;
        subjects = [487:492];
        reps = 1;
    end

end


ROIPath = [StudyPath 'Masks/native_' atlas '/masked/'];


if isequal(win_roi,'y')
    StudyPathnm = ['corr_' atlas '/'];
    smooth = 'n';
    disp('Not smoothing');
    ROIPath_suff = 'peel_';
elseif isequal(output,'r')
    StudyPathnm = ['corr_' atlas '/'];
    smooth = 'y';
    disp('smoothing');
    ROIPath_suff = '';
else
    StudyPathnm = ['corr_z_' atlas '/'];
    smooth = 'y';
    disp('smoothing');
    ROIPath_suff = '';
end




wb_or_n = 'n';
if wb_or_n == 'y'
    StudyPathnm = ['corr_' atlas '_globalsignal/'];
    disp('HAVE A WHOLE-BRAIN REGRESSOR!!');
else
    disp(' NO whole-brain regressor.');
end


cf_or_n = 's';
disp('  Artificially setting use Mia''s correction factor of ~six.');
while cf_or_n ~= 's' && cf_or_n ~= 'c'
    cf_or_n = input('Use Mia''s correction factor of ~six, or calculate (s/c) ?  ','s');
end

block_select = input('Do all blocks? (y/n): ', 's');
if isequal(block_select,'n')
    blocklist = input('Which block(s)? ("1"):  ','s');
else
    blocklist = '123456';
end


count=0;
done = 'n';
while done ~= 'y'
    fnames = input('Please enter the subject number(s) (eg: {"301A"}): ');
    for file=1:length(fnames)
        for l=1:length(subjects)
            if subjects(l) == str2num(fnames{file})
                count=count+1;
                named(count).name = fnames{file};
                %named(count).fnum = fnames{file};
            end
        end
    end
    
    done = input('     Is the number correct (y/n) ?  ','s');
end


for l=1:length(named)
    count2=1;
    for roig = 1:length(ROIgroups)
        for hemi= 1:length(hemis)
            roigname = [ROIgroups{roig} hemis{hemi}];
            
            named(l).group(count2).name = roigname;
            r=dir([ROIPath, ROIPath_suff ROIgroups{roig} '*' hemis{hemi} '*',named(l).name '*.nii']);
        
            count=1;
            for l2=1:length(r)
                if strfind(r(l2).name, named(l).name)
                    thisv = spm_vol([ROIPath, r(l2).name]);
                    [thisY, thisXYZ] = spm_read_vols(thisv);
                    %named(l).roi(count).name = r(l2).name;
                    %named(l).roi(count).idx = find(thisY);
                    named(l).group(count2).roi(count).name = r(l2).name;
                    named(l).group(count2).roi(count).idx =find(thisY);
                    count = count + 1;
                    clear thisY thisXYZ;
                else
                    disp('did not load any rois');
                end
            end
            count2 = count2+1;
        end
    end
end

%load nuisance rois
if regress_nuis
    for l=1:length(named)
        r=dir([NuisPath, '*',named(l).name '*.nii']);
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

clear d r count fnames;     % Remove some unused variables

mask=1;

for l = 1:length(named) %cycle through the subjects
    
    % Load mask to limit voxel analysis space
    disp(' ');
    d = dir([EPIDIR, named(l).name, EPIDIR_suff, named(l).name '-Mean.nii']);
    
    if ~isempty(d)
        disp(['  Loading mean epi ', d(1).name]);
        disp(['    as a simple mask with ', num2str(100*maskthresh),'% threshold...']);
        
        Vmask = spm_vol([EPIDIR,named(l).name, EPIDIR_suff, d(1).name]);
    else
        disp(['  Unable to load mean epi for subject ', named(l).name, '.']);
    end
    
    % Exclude voxels outside the brain
    [thisY, thisXYZ] = spm_read_vols(Vmask);
    if mask
        mthresh = maskthresh * max(max(max(thisY)));
        thisY(find(thisY <= mthresh)) = 0;
        named(l).maskvox = find(thisY);
        named(l).Ysize = size(thisY);
    else
        disp('Not masking against Mean EPI');
        named(l).maskvox = find(thisY);
        named(l).Ysize = size(thisY);
    end
    
    clear thisY thisXYZ;
    
    % Load filenames for EPIS, as well as the motion parameters
    disp(['  Loading epi files for subject ', named(l).name, '...']);
    
    epid=[];
    for e = 1:length(blocklist)
        epid = [epid; dir([EPIDIR,named(l).name,EPIDIR_suff,named(l).name '*' epi_name_1 '*' blocklist(e) '*' epi_name_2])];
    end
    
    cd([EPIDIR,named(l).name,EPIDIR_suff]);
    %%%%%%%%%%%this for loop is new--do each block separate
    for epi=1:length(epid) %cycle through the epis for this subj
        blocklength = length(spm_vol(epid(epi).name))
        actual_epi = blocklist(epi);
        
        named(l).datafiles = [char(ones(length(d(epi)),1) * ...
            [EPIDIR,named(l).name, EPIDIR_suff]),cat(1,epid(epi).name)];
        
        if regress_nuis
            d = dir([EPIDIR,named(l).name, EPIDIR_suff, motionfilename]);
            named(l).motion = load([EPIDIR,named(l).name, EPIDIR_suff, d(epi).name]);
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
        disp('converting rois');
        for l2 = 1:length(named(l).group)
            for l3 = 1:length(named(l).group(l2).roi)
                found = 0;
                notfound = 0;
                named(l).group(l2).roi(l3).mask_space = [];
                named(l).group(l2).roi(l3).mask_missing = [];
                for l4 = 1:length(named(l).group(l2).roi(l3).idx)
                    tmpidx = find(named(l).maskvox == named(l).group(l2).roi(l3).idx(l4));
                    
                    %[I J] = ind2sub(named(l).maskvox == named(l).roi(l3).idx(l4));
                    if ~isempty(tmpidx)
                        found = found+1;
                        named(l).group(l2).roi(l3).mask_space(found) = tmpidx;
                    else
                        notfound = notfound+1;
                        named(l).group(l2).roi(l3).mask_missing(notfound) = named(l).group(l2).roi(l3).idx(l4);
                    end
                end
                if isempty(named(l).group(l2).roi(l3).mask_space)
                    disp(['     --> No voxels surpass threshold for ROI ', named(l).group(l2).roi(l3).name]);
                end
                if ~isempty(named(l).group(l2).roi(l3).mask_missing)
                    disp(['     --> Missing voxels [', int2str(named(l).group(l2).roi(l3).mask_missing), ...
                        '] for ROI ', named(l).group(l2).roi(l3).name]);
                    if length(named(l).group(l2).roi(l3).mask_missing)/length(named(l).group(l2).roi(l3).mask_space) >= 0.75
                        named(l).group(l2).roi(l3).mask_make = 0; %more than 75% of the ROI is outside the EPI
                        disp(['ROI ',named(l).group(l2).roi(l3).name, ' missing more than 75% of voxels!!!!!!!!']);
                    else
                        named(l).group(l2).roi(l3).mask_make = 1;
                    end
                else
                    named(l).group(l2).roi(l3).mask_make = 1;
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
        
        % Apply bandpass filter
        ALL(l).data=detrend(ALL(l).data);
        
        disp('  Applying bandpass filter...');
        for l3 = l:size(ALL(l).data,2)
            ALL(l).data(:,l3) = fbandpass_new(ALL(l).data(:,l3), TR, f, 0);
        end
        
        if isequal(smooth,'y')
            % Smooth filtered data
            disp(['  Smoothing data with Gaussian of FWHM = ', int2str(p_smooth), '...']);
            if ~exist([EPIDIR,named(l).name,EPIDIR_suff,'tmp_epis/'],'dir')
                unix(['mkdir ', EPIDIR,named(l).name, EPIDIR_suff, 'tmp_epis/']);
            end
            
            tmp_epi_dir = [EPIDIR,named(l).name, EPIDIR_suff, 'tmp_epis/'];
            results = Vdata(1);
            for l3 = 1:size(ALL(l).data,1)
                if l3 < 10
                    fnumber = ['00', int2str(l3)];
                elseif l3 < 100
                    fnumber = ['0', int2str(l3)];
                else
                    fnumber = int2str(l3);
                end
                results.fname = [tmp_epi_dir, 'tmp_epis', fnumber, '.nii'];%'.img'];
                xxx = zeros(named(l).Ysize);
                xxx(named(l).maskvox) = ALL(l).data(l3,:);
                
                warning off all;
                spm_write_vol(results, xxx);
                
            end
            for l3 = 1:size(ALL(l).data,1)
                if l3 < 10
                    fnumber = ['00', int2str(l3)];
                elseif l3 < 100
                    fnumber = ['0', int2str(l3)];
                else
                    fnumber = int2str(l3);
                end
                
                warning off all;
                spm_smooth([tmp_epi_dir, 'tmp_epis', fnumber, '.nii'],...
                    [tmp_epi_dir, 'stmp_epis', fnumber, '.nii'],6);
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
        end
        
        % Get whole brain signal
        if wb_or_n == 'y'
            disp('  Computing whole brain signal -- this MUST be better specified...');
            named(l).wholebrain = mean(ALL(l).data,2);
        else
            disp('  NOT including a whole brain signal as a nuisance regressor.');
            named(l).wholebrain = [];
        end
        
        if regress_nuis
            % Get the mean signal for the nuisance variables
            for l3 = 1:length(named(l).nuisroi)
                named(l).nuisroi(l3).data = mean(ALL(l).data(:,named(l).nuisroi(l3).mask_space),2);
            end
        end
        
        % Begin regression to remove components related to nuisance variables
        disp('  Beginning regression of nuisance variables...');
        tlen = size(ALL(l).data,1);
        
        % Diff_nuis approximates the temporal derivative
        
        if regress_nuis
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
        for l2 = 1:length(named(l).group)
            for l3 = 1:length(named(l).group(l2).roi)
                if ~isempty(named(l).group(l2).roi(l3).mask_space)
                    named(l).group(l2).roi(l3).data = mean(ALL(l).data(:,named(l).group(l2).roi(l3).mask_space),2);
                    
                    %this is intended to save only the data within the ROI
                    %(before taking the mean)
                    named(l).group(l2).roi(l3).data2 = ALL(l).data(:,named(l).group(l2).roi(l3).mask_space);
                    
                    if cf_or_n == 'c'
                        named(l).group(l2).roi(l3).correction = sum(xcorr( (named(l).group(l2).roi(l3).data - ...
                            mean(named(l).group(l2).roi(l3).data)), 'coeff').^2);
                    else
                        named(l).group(l2).roi(l3).correction = mias_cf;
                    end
                else
                    named(l).group(l2).roi(l3).data = zeros(size(mean(ALL(l).data(:,named(l).group(l2).roi(l3).mask_space),2)));
                end
            end
        end
        
        
        
        %%%%%%START RUNNING CORRELATIONS%%%%%%
        
        for l5=1:reps
            if reps==2
                blocknum=str2num(actual_epi)*2 - mod(l5,2)
                blocklength=floor(blocklength/2);
            else
                blocknum=str2num(actual_epi);
            end
            
            if l5==1
                start=1;
                stop=blocklength;%Sess(blocknum).dur;
                
            else
                start=blocklength+1;%Sess(blocknum).dur+1;
                stop=start+blocklength-1;%start+Sess(blocknum).dur-1;
            end
            
            if blocknum < 10
                disp('blocknum less');
                epitype = ['Block0' num2str(blocknum)];
            else
                disp('blocknum greater');
                epitype = ['Block' num2str(blocknum)];
            end
            if isequal(vox2vox,'n') || isequal(vox2vox,'b')
                ALL_old = ALL; %save this so do not overwrite
                named_old = named; %save this so do not overwrite
                % Compute the correlations for each of the ROIs with all the other
                % ROIs
                disp('  Computing correlations across all ROIs...');
                
                data=[];
                if blocknum<10
                    strblock=['0' num2str(blocknum)];
                else
                    strblock=num2str(blocknum);
                end
                
                if isequal(win_roi,'y')
                    %within ROI correlations
                    for l2 = 1:length(named(l).group)
                        for l3=1:length(named(l).group(l2).roi)
                            data=[];
                            temproidir = [StudyPath, StudyPathnm, 'corr_'  named(l).group(l2).roi(l3).name(6:end-8) '_mask/'];
                            if ~exist(temproidir,'dir')
                                unix(['mkdir ', temproidir]);
                            end
                            tempfilename = [temproidir, named(l).name, '_Block' strblock];
                            if ~exist([tempfilename '.mat'],'file')
                                disp(['Working on: ', tempfilename ]);
                                for l4=1:size(named(l).roi(l3).mask_space,2)
                                    if named(l).roi(l3).mask_make
                                        data = [data ALL(l).data(start:stop,named(l).group(l2).roi(l3).mask_space(l4))];
                                    else
                                        data = [data zeros(size(ALL(l).data(start:stop,named(l).group(l2).roi(l3).mask_space(l4))))];
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
                    end
                else
                    tempfilename = [StudyPath, StudyPathnm, 'roi2roi/', named(l).name, '_Block' strblock '.mat'];
                    if ~exist(tempfilename)
                        for l2 = 1:length(named(l).group)
                            for l3=1:length(named(l).group(l2).roi)
                                if named(l).group(l2).roi(l3).mask_make
                                    data=[data named(l).group(l2).roi(l3).data(start:stop)];
                                else
                                    data=[data zeros(size(named(l).group(l2).roi(l3).data(start:stop)))];
                                end
                            end
                        end
                        
                        [r] = corr(data);
                        fdata = r;
                        
                        save(tempfilename, 'fdata');
                        clear r fdata;
                    else
                        disp(['exists: ', tempfilename]);
                    end
                    if isequal(partcorr,1)
                        %regress out the contribution of a single ROI
                        for l2 = 1:length(named(l).group)
                            for l4=1:length(named(l).group(l2).roi)
                                partdata = named(l).group(l2).roi(l4).data(start:stop);
                                data=[];
                                for l3=1:length(named(l).group(l2).roi)
                                    if named(l).group(l2).roi(l3).mask_make
                                        data=[data named(l).group(l2).roi(l3).data(start:stop)];
                                    else
                                        data=[data zeros(size(named(l).group(l2).roi(l3).data(start:stop)))];
                                    end
                                end
                                
                                [r] = partialcorr(data,partdata);
                                fdata = r;
                                tempfilename = [StudyPath, StudyPathnm, 'roi2roi/part_', named(l).group(l2).roi(l4).name(1:end-8),'_',named(l).name, '_Block' strblock];
                                save(tempfilename, 'fdata');
                                clear r fdata;
                            end
                        end
                    elseif isequal(partcorr,2)
                        
                        %%%This will regress out the contribution of ALL
                        %%%rois except for the pair
                        for l2 = 1:length(named(l).group)
                            tempfilename = [StudyPath, StudyPathnm, 'roi2roi/partall',named(l).group(l2).name, '_' named(l).name '_Block' strblock '.mat'];
                            
                            if ~exist(tempfilename,'file')
                                %regress out all rois
                                outmat = zeros(length(named(l).group(l2).roi),length(named(l).group(l2).roi));
                                
                                % make an array with all ROI time series
                                alldatavals = [];
                                for l3=1:length(named(l).group(l2).roi)
                                    alldatavals = [alldatavals named(l).group(l2).roi(l3).data(start:stop)];
                                end
                                allroilist = [1:length(named(l).group(l2).roi)];
                                for l3=1:length(named(l).group(l2).roi)
                                    for l4=1:length(named(l).group(l2).roi)
                                        npair1= allroilist(allroilist~=l3);
                                        npair2= allroilist(allroilist~=l4);
                                        nonpairind = intersect(npair1,npair2);
                                        
                                        partdata = alldatavals(:,nonpairind);
                                        
                                        if named(l).group(l2).roi(l3).mask_make && named(l).group(l2).roi(l4).mask_make
                                            data=[named(l).group(l2).roi(l3).data(start:stop) named(l).group(l2).roi(l4).data(start:stop)];
                                        else
                                            data=[zeros(size(named(l).group(l2).roi(l3).data(start:stop))) zeros(size(named(l).group(l2).roi(l3).data(start:stop)))];
                                        end
                                        
                                        outval = partialcorr(data,partdata);
                                        outmat(l3,l4) = outval(1,2);
                                    end
                                end
                                
                                fdata = outmat;
                                save(tempfilename, 'fdata');
                                clear outmat fdata;
                                
                            else
                                disp([tempfilename ' exists']);
                            end
                        end
                   
                    elseif isequal(partcorr,3)
                        
                        %%%This will regress out the contribution of ROIs
                        %%%between 2 streams
                        for l1=1:length(named(l).group)
                            for l2 = 1:length(named(l).group)
                                %if l1 >= l2 || l1==l2
                                %    disp('skip');
                                %else
                                tempfilename = [StudyPath, StudyPathnm, 'roi2roi/partlevel_',named(l).group(l1).name, '_' named(l).group(l2).name, '_' named(l).name '_Block' strblock '.mat'];
                                if ~exist(tempfilename,'file') %~exist(tempfilename,'file')
                                    %regress out all rois
                                    outmat = zeros(length(named(l).group(l1).roi),length(named(l).group(l2).roi));
                                    
                                    % make an array with all ROI time series
                                    alldatavals1 = [];
                                    for l3=1:length(named(l).group(l1).roi)
                                        alldatavals1 = [alldatavals1 named(l).group(l1).roi(l3).data(start:stop)];
                                    end
                                    alldatavals2 = [];
                                    for l3=1:length(named(l).group(l2).roi)
                                        alldatavals2 = [alldatavals2 named(l).group(l2).roi(l3).data(start:stop)];
                                    end
                                    
                                    allroilist1 = [1:length(named(l).group(l1).roi)];
                                    allroilist2 = [1:length(named(l).group(l2).roi)];
                                    
                                    for l3=1:length(named(l).group(l1).roi)
                                        for l4=1:length(named(l).group(l2).roi)
                                            npair1= allroilist1(allroilist1~=l3); %DL1 DL2 DL3
                                            npair2= allroilist2(allroilist2~=l4); %DM1 DM2 DM3
                                            if l1==l2
                                                %within the same stream
                                                %nonpairind = [npair1;npair2];
                                                nonpairind = intersect(npair1,npair2);
                                                partdata = [alldatavals1(:,nonpairind)];
                                            else
                                                partdata = [alldatavals1(:,npair1) alldatavals2(:,npair2)];
                                            end
                                            
                                            if named(l).group(l1).roi(l3).mask_make && named(l).group(l2).roi(l4).mask_make
                                                data=[named(l).group(l1).roi(l3).data(start:stop) named(l).group(l2).roi(l4).data(start:stop)];
                                            else
                                                data=[zeros(size(named(l).group(l1).roi(l3).data(start:stop))) zeros(size(named(l).group(l2).roi(l3).data(start:stop)))];
                                            end
                                            
                                            outval = partialcorr(data,partdata);
                                            outmat(l3,l4) = outval(1,2);
                                        end
                                    end
                                    
                                    fdata = outmat;
                                    %save(tempfilename, 'fdata');
                                    clear outmat fdata;
                                    
                                else
                                    disp([tempfilename ' exists']);
                                end
                            end
                        end
                    elseif isequal(partcorr,4)
                        
                        %%%This will regress out the contribution of ROIs
                        %%%from a single stream (always the first one in
                        %%%the name)
                        for l1=1:length(named(l).group)
                            for l2 = 1:length(named(l).group)
                                %if l1 >= l2 || l1==l2
                                %    disp('skip');
                                %else
                                tempfilename = [StudyPath, StudyPathnm, 'roi2roi/partwithinstream_',named(l).group(l1).name, '_' named(l).group(l2).name, '_' named(l).name '_Block' strblock '.mat'];
                                if ~exist(tempfilename,'file') %~exist(tempfilename,'file')
                                    %regress out all rois
                                    outmat = zeros(length(named(l).group(l1).roi),length(named(l).group(l2).roi));
                                    
                                    % make an array with all ROI time series
                                    alldatavals1 = [];
                                    for l3=1:length(named(l).group(l1).roi)
                                        alldatavals1 = [alldatavals1 named(l).group(l1).roi(l3).data(start:stop)];
                                    end
                                    
                                    allroilist1 = [1:length(named(l).group(l1).roi)];
                                    allroilist2 = [1:length(named(l).group(l2).roi)];
                                    
                                    for l3=1:length(named(l).group(l1).roi)
                                        for l4=1:length(named(l).group(l2).roi)
                                            npair1= allroilist1(allroilist1~=l3); %DL1 DL2 DL3
                                            npair2= allroilist2(allroilist2~=l4); %DM1 DM2 DM3
                                            if l1==l2
                                                %within the same stream
                                                %nonpairind = [npair1;npair2];
                                                nonpairind = intersect(npair1,npair2);
                                                partdata = [alldatavals1(:,nonpairind)];
                                            else
                                                partdata = [alldatavals1(:,npair1)];
                                            end
                                            
                                            if named(l).group(l1).roi(l3).mask_make && named(l).group(l2).roi(l4).mask_make
                                                data=[named(l).group(l1).roi(l3).data(start:stop) named(l).group(l2).roi(l4).data(start:stop)];
                                            else
                                                data=[zeros(size(named(l).group(l1).roi(l3).data(start:stop))) zeros(size(named(l).group(l2).roi(l3).data(start:stop)))];
                                            end
                                            
                                            outval = partialcorr(data,partdata);
                                            outmat(l3,l4) = outval(1,2);
                                        end
                                    end
                                    
                                    fdata = outmat;
                                    save(tempfilename, 'fdata');
                                    clear outmat fdata;
                                    
                                else
                                    disp([tempfilename ' exists']);
                                end
                            end
                        end
                    elseif isequal(partcorr,5)
                        %%%This will regress out the contribution of ROIs
                        %%%in the 'other' stream (not included in the
                        %%%current pair)
                        for l1=1:length(named(l).group)
                            for l2 = 1:length(named(l).group)
                                if l1~= l2
                                    
                                    if l1 == 1 && l2 == 2 || l1 == 2 && l2 == 1
                                        offstream = 3;
                                    elseif l1 == 1 && l2 == 3 || l1 == 3 && l2 == 1
                                        offstream = 2;
                                    elseif l1 == 2 && l2 == 3 || l1 == 3 && l2 == 2
                                        offstream = 1;
                                    end
                                    %if l1 >= l2 || l1==l2
                                    %    disp('skip');
                                    %else
                                    tempfilename = [StudyPath, StudyPathnm, 'roi2roi/partacrossstream_',named(l).group(l1).name, '_' named(l).group(l2).name, '_' named(l).name '_Block' strblock '.mat'];
                                    if ~exist(tempfilename,'file')
                                        %regress out all rois
                                        outmat = zeros(length(named(l).group(l1).roi),length(named(l).group(l2).roi));
                                        
                                        % make an array with all ROI time series
                                        partdata = [];
                                        for l3=1:length(named(l).group(l1).roi)
                                            partdata = [partdata named(l).group(offstream).roi(l3).data(start:stop)];
                                        end
                                        
                                        for l3=1:length(named(l).group(l1).roi)
                                            for l4=1:length(named(l).group(l2).roi)
                                                
                                                if named(l).group(l1).roi(l3).mask_make && named(l).group(l2).roi(l4).mask_make
                                                    data=[named(l).group(l1).roi(l3).data(start:stop) named(l).group(l2).roi(l4).data(start:stop)];
                                                else
                                                    data=[zeros(size(named(l).group(l1).roi(l3).data(start:stop))) zeros(size(named(l).group(l2).roi(l3).data(start:stop)))];
                                                end
                                                
                                                outval = partialcorr(data,partdata);
                                                outmat(l3,l4) = outval(1,2);
                                            end
                                        end
                                        
                                        fdata = outmat;
                                        save(tempfilename, 'fdata');
                                        clear outmat fdata;
                                        
                                    else
                                        disp([tempfilename ' exists']);
                                    end
                                else
                                    disp('l1 = l2');
                                end
                            end
                        end
                    end
                end
            end
            
            if isequal(vox2vox,'y') || isequal(vox2vox,'b')
                
                ALL_old = ALL; %save this so do not overwrite
                named_old = named; %save this so do not overwrite
                
                % Compute the correlations for each of the ROIs across all
                % voxels within the ROIs
                disp('  Computing correlations across all voxels...');
                
                if blocknum<10
                    strblock=['0' num2str(blocknum)];
                else
                    strblock=num2str(blocknum);
                end
                %computing on the whole brain
                for l2 = 1:length(named(l).group)
                    for l3 = 1:length(named(l).group(l2).roi)
                        if ~isempty(named(l).group(l2).roi(l3).mask_space)
                            fileoutname = [StudyPath, StudyPathnm, 'roi2vox/WB_', ...
                                epitype, '_', named(l).group(l2).roi(l3).name];
                            if ~exist(fileoutname,'file')
                                disp(['making: ',fileoutname]);
                                for l4 = 1:size(ALL(l).data,2)
                                    r(l4) = corr(named(l).group(l2).roi(l3).data(start:stop),ALL(l).data(start:stop,l4));
                                end
                                
                                
                                zscore = 0.5 * log ((1+r)./(1-r));
                                zrho = (1/sqrt((tlen-3)/named(l).group(l2).roi(l3).correction));
                                Z = zscore./zrho;
                                
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
                            else
                                disp(['Already made: ', fileoutname]);
                            end
                        end
                    end
                end
            end
            
            ALL = ALL_old; %save this so do not overwrite
            named = named_old; %save this so do not overwrite
            
        end
        
    end
    clear ALL;
end
