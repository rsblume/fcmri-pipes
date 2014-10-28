%function []=make_txt_files(atlas,project_type)

atlas='huettel_all2';
project_type='DFG';
Topdir=['/r/d2/despo/rstate/data/Rest.PFCnet/Data/'];
Corrdir=[Topdir 'corr_' atlas '/roi2roi/'];
cd(Corrdir);

if isequal(project_type,'A')
    subjects = [401:414];
elseif isequal(project_type,'B')
    subjects = [415:420,422:438];
elseif isequal(project_type,'C')
    subjects = [439:445,448:459];
elseif isequal(project_type,'D')
    subjects = [460:465,493:499];
elseif isequal(project_type,'E')
    subjects = [466:484];
elseif isequal(project_type,'F')
    subjects = [485];
elseif isequal(project_type,'G')
    subjects = [487:492];
elseif isequal(project_type,'DFG')
    subjects = [460:465,485,487:499];
end


roigroups = {'DM','FDL','FVL'};
seedval_groups = {[1:6],[1:6],[1:6]};
blocks = {'1'};

%indices of the input matrix 
%DLind_L = [1:2:12];%12
%DLind_R = [2:2:12];
DMoldind_L = [13:2:24];%24
%DMoldind_R = [14:2:24];
DMnewind_L = [25:2:36];
%DMnewind_R = [26:2:36];
DRLind_L = [37:2:46];%10
%DRLind_R = [38:2:46];
FDLind_L = [47:2:58];%12
%FDLind_R = [48:2:58];
FVLind_L = [59:2:70];
%FVLind_R = [60:2:70];
VLind_L = [71:2:78];
VLind_R = [72:2:78];

npDMind=[1:6];
npFDLind=[7:12];
npFVLind=[13:18];


%subnum FDL-FVL DM-FVL DM-FDL
bl = 1;
for sub = 1 :length(subjects)
    for seedvals = 1 :size(seedval_groups,2)
        for seed = 1 :length(seedval_groups{seedvals})
            
            for roi = 1 :size(roigroups,2)
                
                %cycles through ROIs that have been partialed out
                if seed<10 && isequal(atlas,'huettel_all2')
                    filename = ['part_' roigroups{roi} '0' num2str(seed) '_L_' num2str(subjects(sub)) '_Block' blocks{bl}];
                else
                    filename = ['part_' roigroups{roi} num2str(seed) '_L_' num2str(subjects(sub)) '_Block' blocks{bl}];
                end
                load(filename);
                fdata = atanh(fdata);
                
       
                    
                    if isequal(roigroups{roi},'DM') %if partialling DM
                        arr = fdata(FDLind_L,FVLind_L);
                    elseif isequal(roigroups{roi},'FDL') %if partialling FDL
                        arr = fdata(DMoldind_L,FVLind_L);
                    else
                        arr = fdata(DMoldind_L,FDLind_L); %if partialling FVL
                    end
                    
  
                
                
                AllArr(bl,sub,seed,roi) = arr(seed,seed);
                

                % Indices for non-partial fdata don't work on np data
                 npfilename = [num2str(subjects(sub)) '_Block0' blocks{bl}];
                 load(npfilename);
                 fdata = atanh(fdata);
                
                if isequal(roigroups{roi},'DM')
                    arr = fdata(FDLind_L,FVLind_L);
                elseif isequal(roigroups{roi},'FDL')
                    arr = fdata(DMoldind_L,FVLind_L);
                else
                    arr = fdata(DMoldind_L,FDLind_L);
                end
                
                
                AllArr2(bl,sub,seed,roi)=arr(seed,seed);
            end
        end
        
    end
    
end

 save([atlas '_' project_type '_part_Block' blocks{bl} '_fisher'],'AllArr');
 save([atlas '_' project_type '_np_Block' blocks{bl} '_fisher'],'AllArr2');
 
for iroi=1:3;
    for iseed=1:5;
    [x,part_ttests(iroi,iseed)]=ttest(AllArr(1,:,iseed,iroi));
    part_means(iroi,iseed)=mean(AllArr(1,:,iseed,iroi));
    part_std(iroi,iseed)=std(AllArr(1,:,iseed,iroi));
    
    [x,np_ttests(iroi,iseed)]=ttest(AllArr2(1,:,iseed,iroi));
    np_means(iroi,iseed)=mean(AllArr2(1,:,iseed,iroi));
    np_std(iroi,iseed)=std(AllArr2(1,:,iseed,iroi));
    
    
    
    end
end

 
% 
% save([atlas '_' project_type '_non-partial_Block' blocks{bl} '_fisher'],'AllArr2');
% 
% h=figure;
% arr = squeeze(AllArr2);
% boxplot([arr(:,:,1) arr(:,:,2) arr(:,:,3)]);
% xlabel(['FDL-FVL(DM)            ','            DM-FVL(FDL)          ','           DM-FDL(FVL)']);
% title(['non-partial correlation: Block' blocks{1}]);
% saveas(h,[project_type 'orig_boxplot'],'pdf');
% 
% h=figure;
% arr = squeeze(AllArr);
% boxplot([arr(:,:,1) arr(:,:,2) arr(:,:,3)]);
% xlabel(['FDL-FVL(DM)            ','            DM-FVL(FDL)          ','           DM-FDL(FVL)']);
% title(['partial correlation: Block' blocks{1}]);
% saveas(h,[project_type 'partial_boxplot'],'pdf');
% 
% out = mean(AllArr,2);
% out = squeeze(out);
% h=figure; bar(out); legend('FDL-FVL(DM)','DM-FVL(FDL)','DM-FDL(FVL)','Location','Best');
% title(['partial correlation: Block' blocks{1}]);
% xlabel('Level in the hierarchy (caudal -> rostral)'); ylabel('fisher transformed mean correlations');
% %ylim([-1,1]);
% saveas(h,[project_type 'partial'],'pdf');
% out2 = mean(AllArr2,2);
% out2 = squeeze(out2);
% h=figure; bar(out2); legend('FDL-FVL(DM)','DM-FVL(FDL)','DM-FDL(FVL)','Location','Best');
% title(['non-partial correlation: Block' blocks{1}]);
% xlabel('Level in the hierarchy (caudal -> rostral)'); ylabel('fisher transformed mean correlations');
% %ylim([-1,1]);
% saveas(h,[project_type 'orig'],'pdf');
% disp('done');
%fclose(fid);