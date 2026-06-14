# Multi-echo ICA (ME-ICA) Processing of fMRI Data version 4.0.0

# Installation

For a simple install on Linux, open a terminal, set an install directory and run the installer

```bash
INSTALLDIR=$HOME
curl -LO https://raw.githubusercontent.com/ME-ICA/me-ica/3.9.7-gpl/install.sh
bash install.sh $INSTALLDIR
```

This installs all Python dependencies to an independent directory. This installer does not require admin rights to run and does not interfere with existing Python environments. The standard AFNI distribution needs to be installed and available on your path. 

# Usage

Start by activating the environment.

```bash
source ~/activate_meica
```

You can then call `meica.py` from any data directory.

**NOTE**
- ***meica.py* must be *called from inside a data directory*, which will also be where results are stored** 

meica.py has many options which you can view using the -h flag. 

For example:

If ME-fMRI data are called: 		rest_e1.nii.gz rest_e2.nii.gz rest_e3.nii.gz, etc. 
Anatomical is:		mprage.nii.gz

```bash
cd /path/to/data
meica.py -d rest1_e1.nii.gz rest1_e2.nii.gz rest1_e3.nii.gz -a mprage.nii --para_echo --prefix sub1_rest
```

This means:

    -d rest_e1.nii.gz rest_e2...   are the 4-D time series datasets (space separated list) from a multi-echo fMRI acqusition
    -a ...   is a "raw" mprage with a skull
    --para_echo  preprocess echo datasets in parallel for speed
   	--prefix sub1_rest   prefix for final functional output datasets, i.e. sub1_rest_....nii.gz

If you just want to analyze an ME-fMRI run without an anatomical, try something like

```bash
cd /path/to/data
meica.py -d task1_e*.nii.gz
```

N.B. 
- Now a simple space separated list of datasets is supported, no quotes needed as was the case in previous versions. 
- Assuming BIDS .json files are in the folder, TEs are read from the BIDS headers, making usage much easier
- If slice timing information is in BIDS headers, specify a `--tpattern bids` option to `meica.py`. Check AFNI documentation of [3dTshift](http://afni.nimh.nih.gov/pub/dist/doc/program_help/3dTshift.html) to see slice timing codes.

See `meica.py -h` for handling other situations such as: no BIDS headers, anatomical with no skull, no anatomical at all, non-linear warp to standard space, etc.

## Output

Assuming use of `-prefix sub1_rest`,

- `sub1_rest_ctab.txt` : Table of component Kappa, Rho, and variance explained values, plus listing of component classifications.
- `sub1_rest_mefl_nat.nii.gz` : Component maps (in units of \delta S) of ALL ICA components.

Examining accepted components in _ctab.txt and component maps in _mefl.nii.gz (using AFNI or FSL) is key to determining component selection effectiveness.

- `sub1_rest_medn_nat.nii.gz` : 'Denoised' BOLD time series after: basic preprocessing, T2* weighted averaging of echoes (i.e. 'optimal combination'), ICA denoising. Use this dataset for task analysis and resting state time series correlation analysis. See [here](https://www.pnas.org/doi/abs/10.1073/pnas.1301725110) for information on degrees of freedom in denoised data.
- `sub1_rest_tsoc_afw.nii.gz` : 'Raw' BOLD time series dataset after: basic preprocessing and T2* weighted averaging of echoes (i.e. 'optimal combination'). 'Standard' denoising or task analyses can be assessed on this dataset (e.g. motion regression, physio correction, scrubbing, blah...) for comparison to ME-ICA denoising.
- `sub1_rest_mefc_afw.nii.gz` : Component maps (in units of \delta S) of accepted BOLD ICA components. Use this dataset for ME-ICR seed-based connectivity analysis.
- `./meica.rest1_e1/` : contains preprocessing intermediate files.

Suffixes indicate warp spaces (see `meica.py -h`):
- `epi` : Original EPI geometry
- `nat` : Native anatomical affine warp
- `afw` : Affine warp e.g. to MNI space
- `nlw` : Non-linear warp to standard space using Qwarp 

# Some Notes

- For more info on T2* weighted anatomical-functional coregistration click [here](https://pmc.ncbi.nlm.nih.gov/articles/PMC6319659/pdf/nihms-1001520.pdf)
- FWHM smoothing is not recommended. tSNR boost is provided by optimal combination of echoes. For better overlap of 'blobs' across subjects, use non-linear standard space normalization instead with `meica.py ... --qwarp`
