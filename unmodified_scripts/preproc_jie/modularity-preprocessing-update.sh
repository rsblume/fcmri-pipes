

##### copy data from orig folder ###################################3
for subject in 002 003 004 005 006 007 008 009 010 011 012 001
do 
cd /cnari/normal_language/NL-EDS-JY/NL-EDS-SUB/GM-subjects/${subject}
mkdir ./old_res
mv *.* ./old_res
done


for subject in 002 003 004 005 006 007 008 009 010 011 012 001
do 
cd /cnari/normal_language/NL-EDS/subjects/NL-EDS-${subject}/measures/imag-measure/Results
for run in run-1 run-2 run-3
do
cp reg_TS_${run}+orig.* /cnari/normal_language/NL-EDS-JY/NL-EDS-SUB/GM-subjects/${subject}
done
done


for subject in 002 003 004 005 006 007 008 009 010 011 012 001
do 
cd /cnari/normal_language/NL-EDS/subjects/NL-EDS-${subject}/measures/imag-measure/Results
for run in 1 2 3
do
cp motion.${run} /cnari/normal_language/NL-EDS-JY/NL-EDS-SUB/GM-subjects/${subject}
done
done


################ remove wm and ventricles ####################
#### get the mask file ###########

for subject in 006 007 008 009 010 011 012 001 #002 003 004 005 
do 
cd /cnari/normal_language/NL-EDS-JY/scripts/modularity/
python MaskMaker.py --identity ${subject} --volume /cnari/normal_language/NL-EDS-JY/NL-EDS-SUB/${subject}/old_res/volume.ave+orig --automask /cnari/normal_language/NL-EDS-JY/NL-EDS-SUB/${subject}/old_res/automask_Alli+orig --makeautobox y --location /cnari/normal_language/NL-EDS-JY/NL-EDS-SUB/GM-subjects/${subject}/old_res/
done

for subject in 009 010 011 012 001 002 003 004 005 006 007 #008
do 
cd /cnari/normal_language/NL-EDS-JY/NL-EDS-SUB/GM-subjects/${subject}/
for run in run-1 run-2 run-3
do
rm ./old_res/mask.WMclip.frac.${subject}+orig.*
rm ./old_res/mask.VENTclip.frac.${subject}+orig.*
rm ./old_res/mask.WMclip.frac.clip.${subject}+orig.*
rm ./old_res/mask.VENTclip.frac.clip.${subject}+orig.*
3dresample -rmode Li -master reg_TS_run-1+orig -inset ./old_res/mask.WMclip.${subject}+orig -prefix ./old_res/mask.WMclip.frac.${subject}
3dresample -rmode Li -master reg_TS_run-1+orig -inset ./old_res/mask.VENTclip.${subject}+orig -prefix ./old_res/mask.VENTclip.frac.${subject}
3dmerge -1clip 0.9 -prefix ./old_res/mask.WMclip.frac.clip.${subject} ./old_res/mask.WMclip.frac.${subject}+orig.
3dmerge -1clip 0.9 -prefix ./old_res/mask.VENTclip.frac.clip.${subject} ./old_res/mask.VENTclip.frac.${subject}+orig.
3dmaskave -mask ./old_res/mask.WMclip.frac.clip.${subject}+orig reg_TS_${run}+orig | awk '{print $1}' >> wm_${subject}_${run}.1D
3dmaskave -mask ./old_res/mask.VENTclip.frac.clip.${subject}+orig reg_TS_${run}+orig | awk '{print $1}' >> vent${subject}_${run}.1D
done
done

########### smooth data using fwhm 6 #############################

for subject in 010 011 012 001 002 003 004 005 006 007 008 009
do 
cd /cnari/normal_language/NL-EDS-JY/NL-EDS-SUB/GM-subjects/${subject}/
for run in 1 2 3
do
3dmerge -1blur_fwhm 6 -doall -prefix reg_TS_run-${run}_bl6 reg_TS_run-${run}+orig.
done
done


#################################### 3dd get cleanTS #####################################


for subject in 010 011 012 001 002 003 004 005 006 007 008 #009
do 
cd /cnari/normal_language/NL-EDS-JY/NL-EDS-SUB/GM-subjects/${subject}/
for run in 1 2 3
do
3dDeconvolve -jobs 2 -input reg_TS_run-${run}+orig. -polort 3 -tshift -num_stimts 8 \
-stim_file 1 "/cnari/normal_language/NL-EDS-JY/NL-EDS-SUB/GM-subjects/${subject}/motion.${run}[0]" -stim_base 1 -stim_label 1 roll \
-stim_file 2 "/cnari/normal_language/NL-EDS-JY/NL-EDS-SUB/GM-subjects/${subject}/motion.${run}[1]" -stim_base 2 -stim_label 2 pitch \
-stim_file 3 "/cnari/normal_language/NL-EDS-JY/NL-EDS-SUB/GM-subjects/${subject}/motion.${run}[2]" -stim_base 3 -stim_label 3 yaw \
-stim_file 4 "/cnari/normal_language/NL-EDS-JY/NL-EDS-SUB/GM-subjects/${subject}/motion.${run}[3]" -stim_base 4 -stim_label 4 dS \
-stim_file 5 "/cnari/normal_language/NL-EDS-JY/NL-EDS-SUB/GM-subjects/${subject}/motion.${run}[4]" -stim_base 5 -stim_label 5 dL \
-stim_file 6 "/cnari/normal_language/NL-EDS-JY/NL-EDS-SUB/GM-subjects/${subject}/motion.${run}[5]" -stim_base 6 -stim_label 6 dP \
-stim_file 7 "wm_${subject}_run-${run}.1D" -stim_base 7 -stim_label 7 WM \
-stim_file 8 "vent${subject}_run-${run}.1D" -stim_base 8 -stim_label 8 VENT \
-GOFORIT 4 \
-errts cleanTS_${subject}_TS.${run} \
-bucket ${subject}_TS.${run}.buck
done
done

#3dDeconvolve changes:
#-stim_file 15 "wm.1D" -stim_base 15 -stim_label 15 wm \
#-stim_file 16 "vent.1D" -stim_base 16 -stim_label 16 vent \cleanTScat.003.TS.3.corr.nothres
#-errts "error file name"\ #<-This is your cleaned timeseries.

##### delete the first 5 trs, and delete all the rest period between stimuli ###################################

for subject in 003 004 005 006 007 008 009 010 011 012 001 002
do 
cd /cnari/normal_language/NL-EDS-JY/NL-EDS-SUB/GM-subjects/${subject}/
for run in TS.1 TS.2 TS.3
do
#3dTcat -prefix cleanTS_${subject}_${run}_taskonly cleanTS_${subject}_${run}+orig[6..45] cleanTS_${subject}_${run}+orig[67..106] cleanTS_${subject}_${run}+orig[128..167] cleanTS_${subject}_${run}+orig[189..228]   
3dTcat -prefix cleanTS_${subject}_${run}_rest cleanTS_${subject}_${run}+orig[46..60] cleanTS_${subject}_${run}+orig[107..121] cleanTS_${subject}_${run}+orig[168..182] cleanTS_${subject}_${run}+orig[229..243]   
done
done

### concat all rest together as a large baseline

for subject in 003 004 005 006 007 008 009 010 011 012 001 002
do 
cd /cnari/normal_language/NL-EDS-JY/NL-EDS-SUB/GM-subjects/${subject}/
3dTcat -prefix cleanTS_${subject}_restall cleanTS_${subject}_TS.1_rest+orig cleanTS_${subject}_TS.2_rest+orig cleanTS_${subject}_TS.3_rest+orig
done


#################################### 3dd get cleanTS, just remove motion #####################################


for subject in 011 012 001 002 003 004 005 006 007 008 009 #010
do 
cd /cnari/normal_language/NL-EDS-JY/NL-EDS-SUB/GM-subjects/${subject}/
for run in 1 2 3
do
3dDeconvolve -jobs 2 -input reg_TS_run-${run}+orig. -polort 3 -tshift -num_stimts 6 \
-stim_file 1 "/cnari/normal_language/NL-EDS-JY/NL-EDS-SUB/GM-subjects/${subject}/motion.${run}[0]" -stim_base 1 -stim_label 1 roll \
-stim_file 2 "/cnari/normal_language/NL-EDS-JY/NL-EDS-SUB/GM-subjects/${subject}/motion.${run}[1]" -stim_base 2 -stim_label 2 pitch \
-stim_file 3 "/cnari/normal_language/NL-EDS-JY/NL-EDS-SUB/GM-subjects/${subject}/motion.${run}[2]" -stim_base 3 -stim_label 3 yaw \
-stim_file 4 "/cnari/normal_language/NL-EDS-JY/NL-EDS-SUB/GM-subjects/${subject}/motion.${run}[3]" -stim_base 4 -stim_label 4 dS \
-stim_file 5 "/cnari/normal_language/NL-EDS-JY/NL-EDS-SUB/GM-subjects/${subject}/motion.${run}[4]" -stim_base 5 -stim_label 5 dL \
-stim_file 6 "/cnari/normal_language/NL-EDS-JY/NL-EDS-SUB/GM-subjects/${subject}/motion.${run}[5]" -stim_base 6 -stim_label 6 dP \
-GOFORIT 4 \
-errts cleanTS_${subject}_TS_mot.${run} \
-bucket ${subject}_TS_mot.${run}.buck
done
done

#3dDeconvolve changes:
#-stim_file 15 "wm.1D" -stim_base 15 -stim_label 15 wm \
#-stim_file 16 "vent.1D" -stim_base 16 -stim_label 16 vent \cleanTScat.003.TS.3.corr.nothres
#-errts "error file name"\ #<-This is your cleaned timeseries.

##### delete the first 5 trs, and delete all the rest period between stimuli ###################################

for subject in 003 #004 005 006 007 008 009 010 011 012 001 002
do 
cd /cnari/normal_language/NL-EDS-JY/NL-EDS-SUB/GM-subjects/${subject}/
for run in 1 2 3
do
3dTcat -prefix cleanTS_${subject}_TS.${run}_mot_taskonly cleanTS_${subject}_TS_mot.${run}+orig[6..45] cleanTS_${subject}_TS_mot.${run}+orig[67..106] cleanTS_${subject}_TS_mot.${run}+orig[128..167] cleanTS_${subject}_TS_mot.${run}+orig[189..228]   
done
done


####################### try raw data ################################
for subject in 003 #004 005 006 007 008 009 010 011 012 001 002
do 
cd /cnari/normal_language/NL-EDS-JY/NL-EDS-SUB/GM-subjects/${subject}/
for run in TS.1 TS.2 TS.3
do
3dTcat -prefix ${run}_taskonly ${run}+orig[6..45] ${run}+orig[67..106] ${run}+orig[128..167] ${run}+orig[189..228]   
done
done

####################### try reg data ################################

for subject in 003 #004 005 006 007 008 009 010 011 012 001 002
do 
cd /cnari/normal_language/NL-EDS-JY/NL-EDS-SUB/GM-subjects/${subject}/
for run in run-1 run-2 run-3
do
3dTcat -prefix ${run}_reg_taskonly reg_TS_${run}+orig[6..45] reg_TS_${run}+orig[67..106] reg_TS_${run}+orig[128..167] reg_TS_${run}+orig[189..228]   
done
done

