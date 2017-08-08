# Set up script run environment
set -e
export OMP_NUM_THREADS=2
export MKL_NUM_THREADS=2
export AFNI_3dDespike_NEW=YES
mkdir -p meica.sub-001_task-rest_run-01_echo-123_bold
cp _meica_sub-001_task-rest_run-01_echo-123_bold.sh meica.sub-001_task-rest_run-01_echo-123_bold/
cd meica.sub-001_task-rest_run-01_echo-123_bold

# Deoblique, unifize, skullstrip, and/or  autobox anatomical, in starting directory (may take a little while)
if [ ! -e /Users/fmri/Desktop/meica/meica/tests/resources/sub-001_T1w_do.nii.gz ]; then 3dWarp -overwrite -prefix /Users/fmri/Desktop/meica/meica/tests/resources/sub-001_T1w_do.nii.gz -deoblique /Users/fmri/Desktop/meica/meica/tests/resources/sub-001_T1w.nii.gz; fi
if [ ! -e /Users/fmri/Desktop/meica/meica/tests/resources/sub-001_T1w_ns.nii.gz ]; then 3dUnifize -overwrite -prefix /Users/fmri/Desktop/meica/meica/tests/resources/sub-001_T1w_u.nii.gz /Users/fmri/Desktop/meica/meica/tests/resources/sub-001_T1w_do.nii.gz; 3dSkullStrip -shrink_fac_bot_lim 0.3 -orig_vol -overwrite -prefix /Users/fmri/Desktop/meica/meica/tests/resources/sub-001_T1w_ns.nii.gz -input /Users/fmri/Desktop/meica/meica/tests/resources/sub-001_T1w_u.nii.gz; 3dAutobox -overwrite -prefix /Users/fmri/Desktop/meica/meica/tests/resources/sub-001_T1w_ns.nii.gz /Users/fmri/Desktop/meica/meica/tests/resources/sub-001_T1w_ns.nii.gz; fi

# Copy in functional datasets, reset NIFTI tags as needed
3dcalc -a /Users/fmri/Desktop/meica/meica/tests/resources/sub-001_task-rest_run-01_echo-1_bold.nii.gz -expr 'a' -prefix ./sub-001_task-rest_run-01_echo-123_bold.nii
nifti_tool -mod_hdr -mod_field sform_code 1 -mod_field qform_code 1 -infiles ./sub-001_task-rest_run-01_echo-123_bold.nii -overwrite
3dcalc -a /Users/fmri/Desktop/meica/meica/tests/resources/sub-001_task-rest_run-01_echo-2_bold.nii.gz -expr 'a' -prefix ./sub-001_task-rest_run-01_echo-123_bold.nii
nifti_tool -mod_hdr -mod_field sform_code 1 -mod_field qform_code 1 -infiles ./sub-001_task-rest_run-01_echo-123_bold.nii -overwrite
3dcalc -a /Users/fmri/Desktop/meica/meica/tests/resources/sub-001_task-rest_run-01_echo-3_bold.nii.gz -expr 'a' -prefix ./sub-001_task-rest_run-01_echo-123_bold.nii
nifti_tool -mod_hdr -mod_field sform_code 1 -mod_field qform_code 1 -infiles ./sub-001_task-rest_run-01_echo-123_bold.nii -overwrite

# Calculate and save motion and obliquity parameters, despiking first if not disabled, and separately save and mask the base volume
3dWarp -verb -card2oblique ./sub-001_task-rest_run-01_echo-1_bold.nii.gz[0] -overwrite -newgrid 1.000000 -prefix ./sub-001_T1w_ob.nii.gz /Users/fmri/Desktop/meica/meica/tests/resources/sub-001_T1w_ns.nii.gz | \grep  -A 4 '# mat44 Obliquity Transformation ::'  > sub-001_task-rest_run-01_echo-123_bold_obla2e_mat.1D
3dDespike -overwrite -prefix ./sub-001_task-rest_run-01_echo-123_bold_vrA.nii.gz ./sub-001_task-rest_run-01_echo-1_bold.nii.gz 
3daxialize -overwrite -prefix ./sub-001_task-rest_run-01_echo-123_bold_vrA.nii.gz ./sub-001_task-rest_run-01_echo-123_bold_vrA.nii.gz
3dcalc -a ./sub-001_task-rest_run-01_echo-123_bold_vrA.nii.gz[4]  -expr 'a' -prefix eBbase.nii.gz 
3dvolreg -overwrite -tshift -quintic  -prefix ./sub-001_task-rest_run-01_echo-123_bold_vrA.nii.gz -base eBbase.nii.gz -dfile ./sub-001_task-rest_run-01_echo-123_bold_vrA.1D -1Dmatrix_save ./sub-001_task-rest_run-01_echo-123_bold_vrmat.aff12.1D ./sub-001_task-rest_run-01_echo-123_bold_vrA.nii.gz
1dcat './sub-001_task-rest_run-01_echo-123_bold_vrA.1D[1..6]{4..$}' > motion.1D 

# Preliminary preprocessing of functional datasets: despike, tshift, deoblique, and/or axialize

# Preliminary preprocessing dataset sub-001_task-rest_run-01_echo-1_bold.nii.gz of echo 1 to produce e1_ts+orig
3dDespike -overwrite -prefix ./sub-001_task-rest_run-01_echo-1_bold_pt.nii.gz sub-001_task-rest_run-01_echo-1_bold.nii.gz
3dTshift -heptic  -prefix ./e1_ts+orig ./sub-001_task-rest_run-01_echo-1_bold_pt.nii.gz
3drefit -view orig e1_ts*HEAD
3daxialize  -overwrite -prefix ./e1_ts+orig ./e1_ts+orig
3drefit -deoblique -TR 3.0 e1_ts+orig

# Preliminary preprocessing dataset sub-001_task-rest_run-01_echo-2_bold.nii.gz of echo 2 to produce e2_ts+orig
3dDespike -overwrite -prefix ./sub-001_task-rest_run-01_echo-2_bold_pt.nii.gz sub-001_task-rest_run-01_echo-2_bold.nii.gz
3dTshift -heptic  -prefix ./e2_ts+orig ./sub-001_task-rest_run-01_echo-2_bold_pt.nii.gz
3drefit -view orig e2_ts*HEAD
3daxialize  -overwrite -prefix ./e2_ts+orig ./e2_ts+orig
3drefit -deoblique -TR 3.0 e2_ts+orig

# Preliminary preprocessing dataset sub-001_task-rest_run-01_echo-3_bold.nii.gz of echo 3 to produce e3_ts+orig
3dDespike -overwrite -prefix ./sub-001_task-rest_run-01_echo-3_bold_pt.nii.gz sub-001_task-rest_run-01_echo-3_bold.nii.gz
3dTshift -heptic  -prefix ./e3_ts+orig ./sub-001_task-rest_run-01_echo-3_bold_pt.nii.gz
3drefit -view orig e3_ts*HEAD
3daxialize  -overwrite -prefix ./e3_ts+orig ./e3_ts+orig
3drefit -deoblique -TR 3.0 e3_ts+orig

# Prepare T2* and S0 volumes for use in functional masking and (optionally) anatomical-functional coregistration (takes a little while).
3dAllineate -overwrite -final NN -NN -float -1Dmatrix_apply sub-001_task-rest_run-01_echo-3_bold_vrmat.aff12.1D'{4..24}' -base eBbase.nii.gz -input e1_ts+orig'[4..24]' -prefix e1_vrA.nii.gz
3dAllineate -overwrite -final NN -NN -float -1Dmatrix_apply sub-001_task-rest_run-01_echo-3_bold_vrmat.aff12.1D'{4..24}' -base eBbase.nii.gz -input e2_ts+orig'[4..24]' -prefix e2_vrA.nii.gz
3dAllineate -overwrite -final NN -NN -float -1Dmatrix_apply sub-001_task-rest_run-01_echo-3_bold_vrmat.aff12.1D'{4..24}' -base eBbase.nii.gz -input e3_ts+orig'[4..24]' -prefix e3_vrA.nii.gz
3dZcat -prefix basestack.nii.gz e1_vrA.nii.gz e2_vrA.nii.gz e3_vrA.nii.gz
/Users/fmri/miniconda3/bin/python /Users/fmri/Desktop/meica/meica/meica.libs/t2smap.py -d basestack.nii.gz -e 14.5,38.5,62.5
3dUnifize -prefix ./ocv_uni+orig ocv.nii
3dSkullStrip -no_avoid_eyes -prefix ./ocv_ss.nii.gz -overwrite -input ocv_uni+orig
3dcalc -overwrite -a t2svm.nii -b ocv_ss.nii.gz -expr 'a*ispositive(a)*step(b)' -prefix t2svm_ss.nii.gz
3dcalc -overwrite -a s0v.nii -b ocv_ss.nii.gz -expr 'a*ispositive(a)*step(b)' -prefix s0v_ss.nii.gz
3daxialize -overwrite -prefix t2svm_ss.nii.gz t2svm_ss.nii.gz
3daxialize -overwrite -prefix ocv_ss.nii.gz ocv_ss.nii.gz
3daxialize -overwrite -prefix s0v_ss.nii.gz s0v_ss.nii.gz

# Copy anatomical into ME-ICA directory and process warps
cp /Users/fmri/Desktop/meica/meica/tests/resources/sub-001_T1w_ns.nii.gz* .
afnibin_loc=`which 3dSkullStrip`
templateloc=${afnibin_loc%/*}

# If can't find affine-warped anatomical, copy native anatomical here, compute warps (takes a while) and save in start dir. ;  otherwise linkin existing files
if [ ! -e /Users/fmri/Desktop/meica/meica/tests/resources/sub-001_T1w_ns_at.nii.gz ]; then \@auto_tlrc -no_ss -init_xform AUTO_CENTER -base ${templateloc}/MNI_caez_N27+tlrc -input sub-001_T1w_ns.nii.gz -suffix _at
cp sub-001_T1w_ns_at.nii /Users/fmri/Desktop/meica/meica/tests/resources
gzip -f /Users/fmri/Desktop/meica/meica/tests/resources/sub-001_T1w_ns_at.nii
else if [ ! -e /Users/fmri/Desktop/meica/meica/tests/resources/sub-001_T1w_ns_at.nii.gz ]; then ln -s /Users/fmri/Desktop/meica/meica/tests/resources/sub-001_T1w_ns_at.nii.gz .; fi
fi
3dcopy /Users/fmri/Desktop/meica/meica/tests/resources/sub-001_T1w_ns_at.nii.gz sub-001_T1w_ns_at
rm -f sub-001_T1w_ns_at+orig.*; 3drefit -view orig sub-001_T1w_ns_at+tlrc 
3dAutobox -overwrite -prefix ./abtemplate.nii.gz ${templateloc}/MNI_caez_N27+tlrc

# If can't find non-linearly warped anatomical, compute, save back; otherwise link
if [ ! -e /Users/fmri/Desktop/meica/meica/tests/resources/sub-001_T1w_ns_atnl.nii.gz ]; then 

# Compute non-linear warp to standard space using 3dQwarp (get lunch, takes a while) 
3dUnifize -overwrite -GM -prefix ./sub-001_T1w_ns_atu.nii.gz /Users/fmri/Desktop/meica/meica/tests/resources/sub-001_T1w_ns_at.nii.gz
3dQwarp -iwarp -overwrite -resample -useweight -blur 2 2 -duplo -workhard -base ${templateloc}/MNI_caez_N27+tlrc -prefix /Users/fmri/Desktop/meica/meica/tests/resources/sub-001_T1w_ns_atnl.nii.gz -source ./sub-001_T1w_ns_atu.nii.gz
fi
if [ ! -e /Users/fmri/Desktop/meica/meica/tests/resources/sub-001_T1w_ns_atnl.nii.gz ]; then ln -s /Users/fmri/Desktop/meica/meica/tests/resources/sub-001_T1w_ns_atnl.nii.gz .; fi
3daxialize -overwrite -prefix ./sub-001_T1w_ob.nii.gz ./sub-001_T1w_ob.nii.gz
/Users/fmri/miniconda3/bin/python /Users/fmri/Desktop/meica/meica/meica.libs/alignp_mepi_anat.py -t t2svm_ss.nii.gz -a sub-001_T1w_ob.nii.gz -p mepi 
cp alignp.mepi/mepi_al_mat.aff12.1D ./sub-001_T1w_al_mat.aff12.1D
cat_matvec -ONELINE /Users/fmri/Desktop/meica/meica/tests/resources/sub-001_T1w_ns_at.nii.gz::WARP_DATA -I > /Users/fmri/Desktop/meica/meica/tests/resources/sub-001_T1w_xns2at.aff12.1D
cat_matvec -ONELINE /Users/fmri/Desktop/meica/meica/tests/resources/sub-001_T1w_ns_at.nii.gz::WARP_DATA > sub-001_T1w_xat2ns.aff12.1D
cat_matvec -ONELINE sub-001_task-rest_run-01_echo-123_bold_obla2e_mat.1D sub-001_T1w_al_mat.aff12.1D -I > sub-001_task-rest_run-01_echo-123_bold_wmat.aff12.1D
cat_matvec -ONELINE sub-001_task-rest_run-01_echo-123_bold_obla2e_mat.1D sub-001_T1w_al_mat.aff12.1D -I sub-001_task-rest_run-01_echo-123_bold_vrmat.aff12.1D  > sub-001_task-rest_run-01_echo-123_bold_vrwmat.aff12.1D

# Extended preprocessing of functional datasets
3dBrickStat -mask eBbase.nii.gz -percentile 50 1 50 e1_ts+orig[4] > gms.1D
gms=`cat gms.1D`; gmsa=($gms); p50=${gmsa[1]}
voxsize=`ccalc .85*$(3dinfo -voxvol eBbase.nii.gz)**.33`
voxdims="`3dinfo -adi eBbase.nii.gz` `3dinfo -adj eBbase.nii.gz` `3dinfo -adk eBbase.nii.gz`"
echo $voxdims > voxdims.1D
echo $voxsize > voxsize.1D

# Preparing functional masking for this ME-EPI run
3dZeropad  -I 27 -S 27 -A 27 -P 27 -L 27 -R 27  -prefix eBvrmask.nii.gz ocv_ss.nii.gz[0]
3dAllineate -overwrite -final NN -NN -float -1Dmatrix_apply sub-001_task-rest_run-01_echo-123_bold_wmat.aff12.1D -base sub-001_T1w_ns.nii.gz -input eBvrmask.nii.gz -prefix ./eBvrmask.nii.gz -master sub-001_T1w_ns.nii.gz -mast_dxyz 2
3dAllineate -overwrite -final NN -NN -float -1Dmatrix_apply sub-001_task-rest_run-01_echo-123_bold_wmat.aff12.1D -base eBvrmask.nii.gz  -input t2svm_ss.nii.gz -prefix ./t2svm_ss_vr.nii.gz -master sub-001_T1w_ns.nii.gz -mast_dxyz 2
3dAllineate -overwrite -final NN -NN -float -1Dmatrix_apply sub-001_task-rest_run-01_echo-123_bold_wmat.aff12.1D -base eBvrmask.nii.gz -input ocv_uni+orig -prefix ./ocv_uni_vr.nii.gz -master sub-001_T1w_ns.nii.gz -mast_dxyz 2
3dAllineate -overwrite -final NN -NN -float -1Dmatrix_apply sub-001_task-rest_run-01_echo-123_bold_wmat.aff12.1D -base eBvrmask.nii.gz -input s0v_ss.nii.gz -prefix ./s0v_ss_vr.nii.gz -master sub-001_T1w_ns.nii.gz -mast_dxyz 2

# Trim empty space off of mask dataset and/or resample
3dAutobox -overwrite -prefix eBvrmask.nii.gz eBvrmask.nii.gz
3dresample -overwrite -master eBvrmask.nii.gz -dxyz ${voxsize} ${voxsize} ${voxsize} -input eBvrmask.nii.gz -prefix eBvrmask.nii.gz
3dcalc -float -a eBvrmask.nii.gz -expr 'notzero(a)' -overwrite -prefix eBvrmask.nii.gz

# Apply combined co-registration/motion correction parameter set to e1_ts+orig
3dAllineate -final wsinc5 -cubic -float -1Dmatrix_apply sub-001_task-rest_run-01_echo-123_bold_vrwmat.aff12.1D -base eBvrmask.nii.gz -input e1_ts+orig -prefix ./e1_vr.nii.gz
3dTstat -min -prefix ./e1_vr_min.nii.gz ./e1_vr.nii.gz
3dcalc -a eBvrmask.nii.gz -b e1_vr_min.nii.gz -expr 'step(a)*step(b)' -overwrite -prefix eBvrmask.nii.gz
3dcalc -float -overwrite -a eBvrmask.nii.gz -b ./e1_vr.nii.gz[4..$] -expr 'step(a)*b' -prefix ./e1_sm.nii.gz 
3dcalc -float -overwrite -a ./e1_sm.nii.gz -expr "a*10000/${p50}" -prefix ./e1_sm.nii.gz
3dTstat -prefix ./e1_mean.nii.gz ./e1_sm.nii.gz
mv e1_sm.nii.gz e1_in.nii.gz
3dcalc -float -overwrite -a ./e1_in.nii.gz -b ./e1_mean.nii.gz -expr 'a+b' -prefix ./e1_in.nii.gz
3dTstat -stdev -prefix ./e1_std.nii.gz ./e1_in.nii.gz
rm -f e1_pt.nii.gz e1_vr.nii.gz e1_sm.nii.gz

# Apply combined co-registration/motion correction parameter set to e2_ts+orig
3dAllineate -final wsinc5 -cubic -float -1Dmatrix_apply sub-001_task-rest_run-01_echo-123_bold_vrwmat.aff12.1D -base eBvrmask.nii.gz -input e2_ts+orig -prefix ./e2_vr.nii.gz
3dcalc -float -overwrite -a eBvrmask.nii.gz -b ./e2_vr.nii.gz[4..$] -expr 'step(a)*b' -prefix ./e2_sm.nii.gz 
3dcalc -float -overwrite -a ./e2_sm.nii.gz -expr "a*10000/${p50}" -prefix ./e2_sm.nii.gz
3dTstat -prefix ./e2_mean.nii.gz ./e2_sm.nii.gz
mv e2_sm.nii.gz e2_in.nii.gz
3dcalc -float -overwrite -a ./e2_in.nii.gz -b ./e2_mean.nii.gz -expr 'a+b' -prefix ./e2_in.nii.gz
3dTstat -stdev -prefix ./e2_std.nii.gz ./e2_in.nii.gz
rm -f e2_pt.nii.gz e2_vr.nii.gz e2_sm.nii.gz

# Apply combined co-registration/motion correction parameter set to e3_ts+orig
3dAllineate -final wsinc5 -cubic -float -1Dmatrix_apply sub-001_task-rest_run-01_echo-123_bold_vrwmat.aff12.1D -base eBvrmask.nii.gz -input e3_ts+orig -prefix ./e3_vr.nii.gz
3dcalc -float -overwrite -a eBvrmask.nii.gz -b ./e3_vr.nii.gz[4..$] -expr 'step(a)*b' -prefix ./e3_sm.nii.gz 
3dcalc -float -overwrite -a ./e3_sm.nii.gz -expr "a*10000/${p50}" -prefix ./e3_sm.nii.gz
3dTstat -prefix ./e3_mean.nii.gz ./e3_sm.nii.gz
mv e3_sm.nii.gz e3_in.nii.gz
3dcalc -float -overwrite -a ./e3_in.nii.gz -b ./e3_mean.nii.gz -expr 'a+b' -prefix ./e3_in.nii.gz
3dTstat -stdev -prefix ./e3_std.nii.gz ./e3_in.nii.gz
rm -f e3_pt.nii.gz e3_vr.nii.gz e3_sm.nii.gz
3dZcat -overwrite -prefix zcat_ffd.nii.gz  ./e1_in.nii.gz ./e2_in.nii.gz ./e3_in.nii.gz
3dcalc -float -overwrite -a zcat_ffd.nii.gz[0] -expr 'notzero(a)' -prefix zcat_mask.nii.gz

# Perform TE-dependence analysis (takes a good while)
/Users/fmri/miniconda3/bin/python /Users/fmri/Desktop/meica/meica/tedana.py -e 14.5,38.5,62.5 -d zcat_ffd.nii.gz --sourceTEs=-1 --kdaw=10 --rdaw=1  --initcost=tanh --finalcost=tanh --conv=2.5e-5  
voxdims="`3dinfo -adi eBbase.nii.gz` `3dinfo -adj eBbase.nii.gz` `3dinfo -adk eBbase.nii.gz`"
echo $voxdims > voxdims.1D
3dcalc -float -a TED/ts_OC.nii[0] -overwrite -expr 'notzero(a)' -prefix ./export_mask.nii.gz

# Copying results to start directory
3dNwarpApply -overwrite -nwarp /Users/fmri/Desktop/meica/meica/tests/resources/sub-001_T1w_xns2at.aff12.1D '/Users/fmri/Desktop/meica/meica/tests/resources/sub-001_T1w_ns_atnl_WARP.nii.gz' -master abtemplate.nii.gz -dxyz 2 -source TED/ts_OC.nii -interp wsinc5 -prefix sub-001_task-rest_run-01_echo-123_bold_tsoc_nlw.nii 
3dNwarpApply -overwrite -nwarp /Users/fmri/Desktop/meica/meica/tests/resources/sub-001_T1w_xns2at.aff12.1D '/Users/fmri/Desktop/meica/meica/tests/resources/sub-001_T1w_ns_atnl_WARP.nii.gz' -master abtemplate.nii.gz -dxyz 2 -source export_mask.nii.gz -interp wsinc5 -prefix nlw_export_mask.nii 
nifti_tool -mod_hdr -mod_field sform_code 2 -mod_field qform_code 2 -infiles nlw_export_mask.nii -overwrite
3dNotes -h 'T2* weighted average of ME time series' sub-001_task-rest_run-01_echo-123_bold_tsoc_nlw.nii
nifti_tool -mod_hdr -mod_field sform_code 2 -mod_field qform_code 2 -infiles sub-001_task-rest_run-01_echo-123_bold_tsoc_nlw.nii -overwrite
3dcalc -overwrite -a nlw_export_mask.nii -b sub-001_task-rest_run-01_echo-123_bold_tsoc_nlw.nii -expr 'ispositive(a-.5)*b' -prefix sub-001_task-rest_run-01_echo-123_bold_tsoc_nlw.nii ; gzip -f sub-001_task-rest_run-01_echo-123_bold_tsoc_nlw.nii; mv sub-001_task-rest_run-01_echo-123_bold_tsoc_nlw.nii.gz /Users/fmri/Desktop/meica/meica/tests/resources
3dNwarpApply -overwrite -nwarp /Users/fmri/Desktop/meica/meica/tests/resources/sub-001_T1w_xns2at.aff12.1D '/Users/fmri/Desktop/meica/meica/tests/resources/sub-001_T1w_ns_atnl_WARP.nii.gz' -master abtemplate.nii.gz -dxyz 2 -source TED/dn_ts_OC.nii -interp wsinc5 -prefix sub-001_task-rest_run-01_echo-123_bold_medn_nlw.nii 
3dNotes -h 'Denoised timeseries (including thermal noise)' sub-001_task-rest_run-01_echo-123_bold_medn_nlw.nii
nifti_tool -mod_hdr -mod_field sform_code 2 -mod_field qform_code 2 -infiles sub-001_task-rest_run-01_echo-123_bold_medn_nlw.nii -overwrite
3dcalc -overwrite -a nlw_export_mask.nii -b sub-001_task-rest_run-01_echo-123_bold_medn_nlw.nii -expr 'ispositive(a-.5)*b' -prefix sub-001_task-rest_run-01_echo-123_bold_medn_nlw.nii ; gzip -f sub-001_task-rest_run-01_echo-123_bold_medn_nlw.nii; mv sub-001_task-rest_run-01_echo-123_bold_medn_nlw.nii.gz /Users/fmri/Desktop/meica/meica/tests/resources
3dNwarpApply -overwrite -nwarp /Users/fmri/Desktop/meica/meica/tests/resources/sub-001_T1w_xns2at.aff12.1D '/Users/fmri/Desktop/meica/meica/tests/resources/sub-001_T1w_ns_atnl_WARP.nii.gz' -master abtemplate.nii.gz -dxyz 2 -source TED/dn_ts_OC_T1c.nii -interp wsinc5 -prefix sub-001_task-rest_run-01_echo-123_bold_T1c_medn_nlw.nii 
3dNotes -h 'Denoised timeseries with T1 equilibration correction (including thermal noise)' sub-001_task-rest_run-01_echo-123_bold_T1c_medn_nlw.nii
nifti_tool -mod_hdr -mod_field sform_code 2 -mod_field qform_code 2 -infiles sub-001_task-rest_run-01_echo-123_bold_T1c_medn_nlw.nii -overwrite
3dcalc -overwrite -a nlw_export_mask.nii -b sub-001_task-rest_run-01_echo-123_bold_T1c_medn_nlw.nii -expr 'ispositive(a-.5)*b' -prefix sub-001_task-rest_run-01_echo-123_bold_T1c_medn_nlw.nii ; gzip -f sub-001_task-rest_run-01_echo-123_bold_T1c_medn_nlw.nii; mv sub-001_task-rest_run-01_echo-123_bold_T1c_medn_nlw.nii.gz /Users/fmri/Desktop/meica/meica/tests/resources
3dNwarpApply -overwrite -nwarp /Users/fmri/Desktop/meica/meica/tests/resources/sub-001_T1w_xns2at.aff12.1D '/Users/fmri/Desktop/meica/meica/tests/resources/sub-001_T1w_ns_atnl_WARP.nii.gz' -master abtemplate.nii.gz -dxyz 2 -source TED/hik_ts_OC_T1c.nii -interp wsinc5 -prefix sub-001_task-rest_run-01_echo-123_bold_hikts_nlw.nii 
3dNotes -h 'Denoised timeseries with T1 equilibration correction (no thermal noise)' sub-001_task-rest_run-01_echo-123_bold_hikts_nlw.nii
nifti_tool -mod_hdr -mod_field sform_code 2 -mod_field qform_code 2 -infiles sub-001_task-rest_run-01_echo-123_bold_hikts_nlw.nii -overwrite
3dcalc -overwrite -a nlw_export_mask.nii -b sub-001_task-rest_run-01_echo-123_bold_hikts_nlw.nii -expr 'ispositive(a-.5)*b' -prefix sub-001_task-rest_run-01_echo-123_bold_hikts_nlw.nii ; gzip -f sub-001_task-rest_run-01_echo-123_bold_hikts_nlw.nii; mv sub-001_task-rest_run-01_echo-123_bold_hikts_nlw.nii.gz /Users/fmri/Desktop/meica/meica/tests/resources
3dNwarpApply -overwrite -nwarp /Users/fmri/Desktop/meica/meica/tests/resources/sub-001_T1w_xns2at.aff12.1D '/Users/fmri/Desktop/meica/meica/tests/resources/sub-001_T1w_ns_atnl_WARP.nii.gz' -master abtemplate.nii.gz -dxyz 2 -source TED/betas_hik_OC.nii -interp wsinc5 -prefix sub-001_task-rest_run-01_echo-123_bold_mefc_nlw.nii 
3dNotes -h 'Denoised ICA coeff. set for ME-ICR seed-based FC analysis' sub-001_task-rest_run-01_echo-123_bold_mefc_nlw.nii
nifti_tool -mod_hdr -mod_field sform_code 2 -mod_field qform_code 2 -infiles sub-001_task-rest_run-01_echo-123_bold_mefc_nlw.nii -overwrite
3dcalc -overwrite -a nlw_export_mask.nii -b sub-001_task-rest_run-01_echo-123_bold_mefc_nlw.nii -expr 'ispositive(a-.5)*b' -prefix sub-001_task-rest_run-01_echo-123_bold_mefc_nlw.nii ; gzip -f sub-001_task-rest_run-01_echo-123_bold_mefc_nlw.nii; mv sub-001_task-rest_run-01_echo-123_bold_mefc_nlw.nii.gz /Users/fmri/Desktop/meica/meica/tests/resources
3dNwarpApply -overwrite -nwarp /Users/fmri/Desktop/meica/meica/tests/resources/sub-001_T1w_xns2at.aff12.1D '/Users/fmri/Desktop/meica/meica/tests/resources/sub-001_T1w_ns_atnl_WARP.nii.gz' -master abtemplate.nii.gz -dxyz 2 -source TED/betas_OC.nii -interp wsinc5 -prefix sub-001_task-rest_run-01_echo-123_bold_mefl_nlw.nii 
3dNotes -h 'Full ICA coeff. set for component assessment' sub-001_task-rest_run-01_echo-123_bold_mefl_nlw.nii
nifti_tool -mod_hdr -mod_field sform_code 2 -mod_field qform_code 2 -infiles sub-001_task-rest_run-01_echo-123_bold_mefl_nlw.nii -overwrite
3dcalc -overwrite -a nlw_export_mask.nii -b sub-001_task-rest_run-01_echo-123_bold_mefl_nlw.nii -expr 'ispositive(a-.5)*b' -prefix sub-001_task-rest_run-01_echo-123_bold_mefl_nlw.nii ; gzip -f sub-001_task-rest_run-01_echo-123_bold_mefl_nlw.nii; mv sub-001_task-rest_run-01_echo-123_bold_mefl_nlw.nii.gz /Users/fmri/Desktop/meica/meica/tests/resources
3dNwarpApply -overwrite -nwarp /Users/fmri/Desktop/meica/meica/tests/resources/sub-001_T1w_xns2at.aff12.1D '/Users/fmri/Desktop/meica/meica/tests/resources/sub-001_T1w_ns_atnl_WARP.nii.gz' -master abtemplate.nii.gz -dxyz 2 -source TED/feats_OC2.nii -interp wsinc5 -prefix sub-001_task-rest_run-01_echo-123_bold_mefcz_nlw.nii 
3dNotes -h 'Z-normalized spatial component maps' sub-001_task-rest_run-01_echo-123_bold_mefcz_nlw.nii
nifti_tool -mod_hdr -mod_field sform_code 2 -mod_field qform_code 2 -infiles sub-001_task-rest_run-01_echo-123_bold_mefcz_nlw.nii -overwrite
3dcalc -overwrite -a nlw_export_mask.nii -b sub-001_task-rest_run-01_echo-123_bold_mefcz_nlw.nii -expr 'ispositive(a-.5)*b' -prefix sub-001_task-rest_run-01_echo-123_bold_mefcz_nlw.nii ; gzip -f sub-001_task-rest_run-01_echo-123_bold_mefcz_nlw.nii; mv sub-001_task-rest_run-01_echo-123_bold_mefcz_nlw.nii.gz /Users/fmri/Desktop/meica/meica/tests/resources
cp TED/comp_table.txt /Users/fmri/Desktop/meica/meica/tests/resources/sub-001_task-rest_run-01_echo-123_bold_ctab.txt
cp TED/meica_mix.1D /Users/fmri/Desktop/meica/meica/tests/resources/sub-001_task-rest_run-01_echo-123_bold_mmix.1D