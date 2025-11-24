# Multi-echo ICA (ME-ICA) Processing of fMRI Data


# Dependencies

1. AFNI
2. Python 3.9
3. numpy
4. scipy


# Installation

Install Python and other dependencies. If you have AFNI installed and on your path, you should already have an up-to-date version of ME-ICA on your path. Running `meica.py` without any options will check for dependencies and let you know if they are met. If you don't have numpy/scipy (or appropriate versions) installed, I would strongly recommend using the [Enthought Canopy Python Distribution](https://www.enthought.com/downloads/). Click [here](http://wiki.org/installing.html) for more installation help.

# Important Files and Directories

- `meica.py` : a master script that performs preprocessing and calls the ICA/TE-dependence analysis script `lib_tedana.py`
- `libmeica` : a folder that includes utility functions for TE-dependence analysis for denoising and anatomical-functional co-registration
- `tedana.py` : performs ICA and TE-dependence calculations

# Usage

fMRI data is called: 		rest_e1.nii.gz rest_e2.nii.gz rest_e3.nii.gz, etc. 
Anatomical is:		mprage.nii.gz

meica.py and tedana.py have a number of options which you can view using the -h flag. 

Here's an example use:

    meica.py -d rest1_e1.nii.gz rest1_e2.nii.gz rest1_e3.nii.gz  -a mprage.nii --MNI --prefix sub1_rest

This means:

    -d rest_e1.nii.gz rest_e2...   are the 4-D time series datasets (space separated list) from a multi-echo fMRI acqusition
    -a ...   is a "raw" mprage with a skull
    --MNI   warp anatomical to MNI space using a built-in high-resolution MNI template. 
   	--prefix sub1_rest   prefix for final functional output datasets, i.e. sub1_rest_....nii.gz

N.B. 
- Now a simple space separated list of datasets is supported, no quotes needed. 
- Assuming BIDS .json files are in the folder, TEs are read from the BIDS headers, making usage much easier

Other options:
    -e 15,30,45   are the echo times in milliseconds

Again, see `meica.py -h` for handling other situations such as: anatomical with no skull, no anatomical at all, applying FWHM smoothing, non-linear warp to standard space, etc.

Click [here](http://wiki.org/group_analysis.html) more info on group analysis.

## Output

- `./meica.rest1_e1/` : contains preprocessing intermediate files. Click [here](http://wiki.org/meica_preprocessing.html) for detailed listing.
- `sub1_rest_medn.nii.gz` : 'Denoised' BOLD time series after: basic preprocessing, T2* weighted averaging of echoes (i.e. 'optimal combination'), ICA denoising. Use this dataset for task analysis and resting state time series correlation analysis. See [here](http://wiki.org/viewing_results.html#dof) for information on degrees of freedom in denoised data.
- `sub1_rest_tsoc.nii.gz` : 'Raw' BOLD time series dataset after: basic preprocessing and T2* weighted averaging of echoes (i.e. 'optimal combination'). 'Standard' denoising or task analyses can be assessed on this dataset (e.g. motion regression, physio correction, scrubbing, blah...) for comparison to ME-ICA denoising.
- `sub1_rest_mefc.nii.gz` : Component maps (in units of \delta S) of accepted BOLD ICA components. Use this dataset for ME-ICR seed-based connectivity analysis.
- `sub1_rest_mefl.nii.gz` : Component maps (in units of \delta S) of ALL ICA components.
- `sub1_rest_ctab.nii.gz` : Table of component Kappa, Rho, and variance explained values, plus listing of component classifications. See [here](http://wiki.org/viewing_results.html#kappa_spectra) for more info.

For a step-by-step guide on how to assess ME-ICA results in more detail, click [here](http://wiki.org/viewing_results.html)

#Some Notes

- Make sure your datasets have slice timing information in the header. If slice timing information is in BIDS headers, specify a `--tpattern bids` option to `meica.py`. Check AFNI documentation of [3dTshift](http://afni.nimh.nih.gov/pub/dist/doc/program_help/3dTshift.html) to see slice timing codes.
- For more info on T2* weighted anatomical-functional coregistration click [here](http://wiki.org/meica_alignp_mepi_anat.html)
- FWHM smoothing is not recommended. tSNR boost is provided by optimal combination of echoes. For better overlap of 'blobs' across subjects, use non-linear standard space normalization instead with `meica.py ... --qwarp`
