### To use these parcellations, first convert them into BRIK files for AfNI analysis: 

### try sub12 first, use orignial spec file, and 2005 / 2009 annot files

for subject in 001 002 003 004 005 006 007 008 009 010 011 # 012
do
cd /cnari/normal_language/NL-EDS/freesurf/subjects/${subject}/freesurfer/SUMA
for hemi in lh rh
do
3dSurf2Vol \
-spec Freesurfer_${hemi}.spec \
-surf_A smoothwm \
-grid_parent /cnari/normal_language/NL-EDS-JY/NL-EDS-SUB/GM-subjects/${subject}/cleanTS_${subject}_TS.2+orig \
-sv anatomy_noskull_AlndExp+orig \
-map_func max \
-prefix ./${hemi}_parc_raw \
-sdata_1D ${hemi}.aparc.a2005s.annot.1D.roi
done
done


#### interaction effect for Context and Reality ######

## left IPG/AG: 473
## right IPG/AG: 473
## left IPS: 510
## left precuneus: 478

for subject in 0612 0613 0643 0644 0645 0637 0641 0653 0634 0642 0615 17SUB 0649 0651 0638 0639 0650 0640 0648 0627 
do
cd /cnari/normal_language/NL-FTC/FTC${subject}/results/jie_try/freesurfer/SUMA
for num in 473 #473 510 478 
do
for hemi in rh #lh
do
#rm ${hemi}_${num}_anat_parc+orig.*
3dcalc -prefix ${hemi}_${num}_anat_parc -a ${hemi}_parc_raw+orig -expr "step(equals(a,${num}))"
done
done
done


##### ??? the RH roi mask is empty.... OK, problem solved... it has to be the same number with the left hemi... why???

for subject in 0612 0613 0643 0644 0645 0637 0641 0653 0634 0642 0615 17SUB 0649 0651 0638 0639 0650 0640 0648 0627 
do
cd /cnari/normal_language/NL-FTC/FTC${subject}/results/jie_try/freesurfer/SUMA
for num in 473 510 478 # 555
do
for hemi in lh #rh
do
3dmaskave -mask ${hemi}_${num}_anat_parc+orig -quiet ../../${subject}_ftc_context_onset_dur_dmblock_mulp1_3dd_3ddres+orig.[21] > ${subject}.roi_human1st_${hemi}_${num}.txt
3dmaskave -mask ${hemi}_${num}_anat_parc+orig -quiet ../../${subject}_ftc_context_onset_dur_dmblock_mulp1_3dd_3ddres+orig.[24] > ${subject}.roi_human3rd_${hemi}_${num}.txt
3dmaskave -mask ${hemi}_${num}_anat_parc+orig -quiet ../../${subject}_ftc_context_onset_dur_dmblock_mulp1_3dd_3ddres+orig.[27] > ${subject}.roi_unreal1st_${hemi}_${num}.txt
3dmaskave -mask ${hemi}_${num}_anat_parc+orig -quiet ../../${subject}_ftc_context_onset_dur_dmblock_mulp1_3dd_3ddres+orig.[30] > ${subject}.roi_unreal3rd_${hemi}_${num}.txt
done
done
done

