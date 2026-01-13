# Multi-echo ICA (ME-ICA) Processing of fMRI Data version 3.9.7

# Installation

For a simple install on Linux, open a terminal, `cd` to a directory of your choosing (e.g. `$HOME`), and clone this repository. 

```bash
cd $HOME
git clone -b 3.9.7-gpl --single-branch https://github.com/ME-ICA/me-ica me-ica-3.9.7
```

Run the installer, specifying the installation directory (e.g. could be `$HOME` or `~/bin`).

```bash
INSTALLATION_DIR=$HOME
bash $HOME/me-ica-3.9.7/install.sh $INSTALLATION_DIR
```

You must already have AFNI installed and available on the path (e.g. try executing `to3d` on the command line).

# Usage

Start by activating the environment. Assuming your INSTALLATION_DIR was $HOME

```bash
source ~/activate_meica
```

You can then call `meica.py` from any location.

If ME-fMRI data are called: 		rest_e1.nii.gz rest_e2.nii.gz rest_e3.nii.gz, etc. 
Anatomical is:		mprage.nii.gz

meica.py has a number of options which you can view using the -h flag. 

**NOTE**
- ***meica.py* (and this codebase) can be *installed to any directory* or be on the path**
- ***meica.py* must be *called from inside the data directory*, which will also be where results are stored** 

For example:

```bash
cd /path/to/data
meica.py -d rest1_e1.nii.gz rest1_e2.nii.gz rest1_e3.nii.gz  -a mprage.nii --MNI --prefix sub1_rest
```

This means:

    -d rest_e1.nii.gz rest_e2...   are the 4-D time series datasets (space separated list) from a multi-echo fMRI acqusition
    -a ...   is a "raw" mprage with a skull
    --MNI   warp anatomical to MNI space using a built-in high-resolution MNI template. 
   	--prefix sub1_rest   prefix for final functional output datasets, i.e. sub1_rest_....nii.gz

If you just want to analyze an ME-fMRI without an anatomical, do

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

- `./meica.rest1_e1/` : contains preprocessing intermediate files.
- `sub1_rest_medn.nii.gz` : 'Denoised' BOLD time series after: basic preprocessing, T2* weighted averaging of echoes (i.e. 'optimal combination'), ICA denoising. Use this dataset for task analysis and resting state time series correlation analysis. See [here](http://wiki.org/viewing_results.html#dof) for information on degrees of freedom in denoised data.
- `sub1_rest_tsoc.nii.gz` : 'Raw' BOLD time series dataset after: basic preprocessing and T2* weighted averaging of echoes (i.e. 'optimal combination'). 'Standard' denoising or task analyses can be assessed on this dataset (e.g. motion regression, physio correction, scrubbing, blah...) for comparison to ME-ICA denoising.
- `sub1_rest_mefc.nii.gz` : Component maps (in units of \delta S) of accepted BOLD ICA components. Use this dataset for ME-ICR seed-based connectivity analysis.
- `sub1_rest_mefl.nii.gz` : Component maps (in units of \delta S) of ALL ICA components.
- `sub1_rest_ctab.nii.gz` : Table of component Kappa, Rho, and variance explained values, plus listing of component classifications.

# Some Notes

- For more info on T2* weighted anatomical-functional coregistration click [here](https://pmc.ncbi.nlm.nih.gov/articles/PMC6319659/pdf/nihms-1001520.pdf)
- FWHM smoothing is not recommended. tSNR boost is provided by optimal combination of echoes. For better overlap of 'blobs' across subjects, use non-linear standard space normalization instead with `meica.py ... --qwarp`
