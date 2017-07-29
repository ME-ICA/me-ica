# ME-QC

Quality control metrics for data processed using Prantik Kundu's [Multi-Echo Independent Components Analysis (ME-ICA)](https://bitbucket.org/prantikk/me-ica/overview). 

[![Build Status](https://travis-ci.org/emdupre/meqc.svg?branch=master)](https://travis-ci.org/emdupre/meqc) [![Coverage Status](https://coveralls.io/repos/github/emdupre/meqc/badge.svg?branch=master)](https://coveralls.io/github/emdupre/meqc?branch=master)

## Usage

ME-ICA minimally requires (1) acquired echo times (in milliseconds) and (2) functional datasets equal to the number of acquired echoes. But you can supply many other options, viewable with `meica.py -h`.  

Here's one example use case:

    meica.py -d 'sub-001_task-rest_echo-[1,2,3]_run-01_meep1.nii.gz' -e 15,30,45 -b 12s -a 'sub-001_T1w.nii.gz' --MNI

Where:

`-e 15,30,45`  are the echo times in milliseconds  
`-d 'sub-001_task-rest_echo-[1,2,3]_run-01_meep1.nii.gz'`  are the 4-D multi-echo fMRI datasets. Can supply each dataset individually or with bash shell expansion  
`-a 'sub-001_T1w.nii.gz'`  is an anatomical image with skull  
`-b 12s`  means drop first 12 seconds of data for equilibration  
`--MNI`  warp anatomical to MNI space using AFNI's MNI_caez_N27 (Colin 27) template. 

Additional, optional parameters support situations such as: anatomical with no skull, applying FWHM smoothing, non-linear warping, etc.

## Procedure
PROCEDURE 1 : Preprocess multi-echo datasets and apply multi-echo ICA based
on spatial concatenation
-Check arguments, input filenames, and filesystem for dependencies
-Calculation of motion parameters based on images with highest constrast
-Application of motion correction and coregistration parameters
-Misc. EPI preprocessing (temporal alignment, smoothing, etc) in appropriate
order
-Compute PCA and ICA in conjuction with TE-dependence analysis

## Output

- `*_medn.nii.gz` : 'Denoised' BOLD time series after: basic preprocessing, T2* weighted averaging of echoes (i.e. 'optimal combination'), ICA denoising. Use this dataset for task analysis and resting state time series correlation analysis.
- `*_tsoc.nii.gz` : 'Raw' BOLD time series dataset after: basic preprocessing and T2* weighted averaging of echoes (i.e. 'optimal combination'). 'Standard' denoising or task analyses can be assessed on this dataset (e.g. motion regression, physio correction, scrubbing, blah...) for comparison to ME-ICA denoising.
- `*_mefc.nii.gz` : Component maps (in units of \delta S) of accepted BOLD ICA components. Use this dataset for ME-ICR seed-based connectivity analysis.
- `*_mefl.nii.gz` : Component maps (in units of \delta S) of ALL ICA components.
- `*_ctab.nii.gz` : Table of component Kappa, Rho, and variance explained values, plus listing of component classifications. 

## Some Notes

- Make sure your datasets have slice timing information in the header. If not sure, specify a `--tpattern` option to `meica.py`. Check AFNI documentation of [3dTshift](http://afni.nimh.nih.gov/pub/dist/doc/program_help/3dTshift.html) to see slice timing codes.
- FWHM smoothing is not recommended. tSNR boost is provided by optimal combination of echoes. For better overlap of 'blobs' across subjects, use non-linear standard space normalization instead with `meica.py ... --qwarp`

## Citations

If you use ME-ICA in publications, please cite:

> Kundu, P., Brenowitz, N.D., Voon, V., Worbe, Y., Vertes, P.E., Inati, S.J.,
Saad, Z.S., Bandettini, P.A. & Bullmore, E.T. Integrated strategy for
improving functional connectivity mapping using multiecho fMRI. PNAS (2013).
http://dx.doi.org/10.1073/pnas.1301725110

> Kundu, P., Inati, S.J., Evans, J.W., Luh, W.M. & Bandettini, P.A.
Differentiating BOLD and non-BOLD signals in fMRI time series using
multi-echo EPI. NeuroImage (2011).
http://dx.doi.org/10.1016/j.neuroimage.2011.12.028
