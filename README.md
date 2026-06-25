# Multi-Echo ICA (MEICA4) Processing of fMRI Data

## Version 4.0.1

Multi-Echo ICA v4, or MEICA4, is a pipeline for comprehensive processing of multi-echo fMRI data. MEICA4 uses the echo-time dependence of BOLD signal changes to distinguish BOLD-like signal sources from non-BOLD sources. Other mechanism-aware measures are also used to identify and remove complex artifacts. This approach provides robust denoising and analysis of multi-echo fMRI datasets from task, resting-state, and ecologically valid paradigms across a wide range of human-subject and patient-research studies, using a single pipeline with uniform configuration. *Note ME-ICA is a research tool and not for standard clinical use.*

MEICA4 is designed to reduce the practical barrier to entry for multi-echo fMRI. On Linux and Windows (via Windows Subsystem for Linux, WSL2), the standard installer downloads a matched ME-ICA AFNI runtime into the MEICA4 source directory. This runtime supplies the AFNI command-line tools, templates, runtime libraries, and selected utilities needed by MEICA4 without requiring users to install AFNI separately or make system-level changes.

# Installation

For a user installation on Linux or WSL2 without administrator privileges, copy, paste, and run this one-line command in a terminal. This assumes installation into your home directory:

```bash
curl -fsSL https://raw.githubusercontent.com/ME-ICA/me-ica/4.0.1/install.sh | bash -s -- "$HOME/meica4"
```

This downloads and executes `install.sh`, which installs MEICA4 into a dedicated micromamba Python environment that does not interfere with other Python installations. On Linux and WSL2, the installer also downloads a matched ME-ICA AFNI runtime release and installs it inside the MEICA4 source directory.

macOS users must install AFNI separately and ensure AFNI commands such as `3dinfo` and `3dSkullStrip` are available on `PATH` before running the installer.

# Usage

After using the installer, start usage by activating the MEICA4 environment:

```bash
source ~/activate_meica
```

You can then call `meica.py` from any directory.

One approach is to run `meica.py` from inside the data directory. For example, given functional and anatomical NIfTI data that have BIDS `.json` sidecars, for complete preprocessing, denoising, affine warp to standard MNI space, and output of native subject-space equivalents:

```bash
cd /path/to/data
meica.py -d rest1_e*.nii.gz -a mprage.nii.gz --MNI --native
```

You can also run `meica.py` from another directory by giving full or relative paths to the input datasets. For example, given data in BIDS organization, and for nonlinear warp to MNI space, slice timing correction, and all files output with a set prefix:

```bash
cd /path/to/results
meica.py \
  -d func/sub-01_ses-01_task-rest_e1.nii.gz \
     func/sub-01_ses-01_task-rest_e2.nii.gz \
     func/sub-01_ses-01_task-rest_e3.nii.gz \
  -a /path/to/anat/mprage.nii.gz \
  --tpattern bids \
  --MNI \
  --qwarp \
  --prefix sub1_rest
```

In this example, make sure the BIDS sidecars (`.json`) include slice-timing information. The outputs of the run will be stored in the directory where MEICA4 was called from, `/path/to/results` in this case.

On Linux and WSL2 installations that include the distributed runtime, activation also makes `dcm2niix` available. For example, to convert a folder of DICOMs to NIfTI:

```bash
mkdir -p nifti_out
dcm2niix -z y -b y -o nifti_out /path/to/dicom_folder
```


## Minimal use

MEICA4 does not require an anatomical image for analysis. Assuming NIfTI files and BIDS sidecars are available, the simplest use of MEICA4 for comprehensive isolation of BOLD signals from artifacts is:

```bash
meica.py -d /path/to/data/rest1_e*.nii.gz
```

This command produces outputs with an `_epi.nii.gz` suffix (see below).

If BIDS sidecars are not present, provide echo times explicitly with `-e`:

```bash
meica.py -d /path/to/data/rest1_e*.nii.gz -e 10 20 30
```

Co-registration to an anatomical image processed by a different pipeline can then be performed using MEICA4 outputs, especially the optimally combined output image, `*_tsoc_epi.nii.gz`; see below for output details.

## Data inputs

MEICA4 is designed around modern NIfTI and BIDS-style metadata: imaging inputs ending in `.nii` or `.nii.gz`, accompanied by `.json` sidecars that hold key acquisition metadata. For scan parameters to be read from the data, one BIDS `.json` sidecar must be present for each NIfTI dataset, using the same prefix, for example `task1_echo1.nii.gz` and `task1_echo1.json`.

If BIDS `.json` files are not available, echo times must be specified explicitly, for example using `-e 10 20 30` corresponding to echo times of 10, 20, and 30 milliseconds.

For older or archival data, MEICA4 will still accept inputs in AFNI `+orig` format, for example:

```bash
meica.py -d task1-echo*+orig.HEAD -e 10 20 30
```

## Performance

MEICA4 is engineered for high-speed processing and efficient use of system memory through extensive use of parallelism and disk-backed memory mapping. The intention is to support expedient analysis of many-subject datasets from large studies as well as near-line use in clinical research where an intervention is conducted on the basis of results.

Use the option `-j N` to change the number of CPUs made available for analysis. The default is `-j 4`. Use `-j 1` if strict single-CPU behavior is desired. Call as `/usr/bin/time meica.py ...` to measure performance and timing of a run.

On a modern laptop, the examples above can often complete analysis in approximately the time required to acquire the data. For data with spatial resolution around 3 mm and TR around 2 seconds, about 10 minutes of processing time may be needed for functional and anatomical data that take about 20 minutes to acquire.

Analyzing data acquired at around 2 mm isotropic resolution and TR < 2 seconds over approximately 15 minutes can require around 16 GB of memory. For best performance, conduct analysis on solid-state storage.

## Disk usage

ME-fMRI data and processing intermediates can occupy a large amount of disk space. If a sizeable temporary disk or faster write performance on another (i.e. solid state) disk is available, use `--work_root DIR` to conduct processing in another directory and copy final data files back to the start directory when complete.

The default behavior is to delete most intermediate files to save disk space while preserving key intermediates for later inspection. For deeper debugging, use `--keep_int` to retain preprocessing intermediates in the `meica.<run_label>/` processing directory. Use `--no_trim_work` to retain the full working directory produced by the denoising processing.

## Detailed options

```text
-d DATASETS
    Multi-echo fMRI datasets, given as a space-separated list.
    Example: -d rest_e1.nii.gz rest_e2.nii.gz rest_e3.nii.gz

-e TE1 TE2 TE3 ...
    Echo times in milliseconds.
    Default is BIDS metadata lookup.

-a ANAT
    Anatomical dataset. Optional.

-b TIME
    Time to steady-state equilibration, specified in seconds or volumes.
    Examples: -b 10s or -b 5v

--prefix PREFIX
    Prefix for final MEICA4 output datasets.

-j N / --cpus N
    Maximum number of CPUs/OpenMP threads to use.

--MNI
    Warp to MNI standard space using the default high-resolution template.

--space TEMPLATE
    Use a specific standard-space template.

--qwarp
    Use nonlinear warping to standard space.
    Requires an anatomical dataset and a standard-space template, for example with --MNI or --space.

--native
    Output native-space results in addition to standard-space results.

--fres RES
    Specify isotropic functional voxel size for resampling during preprocessing.
    Example: --fres 2.5

--tpattern TPATTERN
    Slice timing pattern.
    Use --tpattern bids to read slice timing from BIDS metadata.

--no_skullstrip
    Use when the anatomical dataset is already skull-stripped and intensity-normalized.

--no_despike
    Disable functional despiking.
    Despiking is enabled by default and recommended.

--smooth FWHM
    Apply spatial smoothing.
    Default is off.

--work_root DIR
    Write the MEICA4 working directory somewhere other than the current directory during processing.

--no_trim_work
    Retain the full working directory after processing.

--keep_int
    Keep preprocessing intermediate files.

--script_only
    Generate the processing script only, then exit.

--pp_only
    Run preprocessing only, then exit.

--test_proc
    Align and preprocess one dataset, then exit.
    Useful for testing.

--skip_check
    Skip dependency checks during initialization.

--RESUME
    Attempt to resume from normalization onward.

--OVERWRITE
    Overwrite an existing MEICA4 output directory after a short warning delay.

--OVERWRITE_NOWAIT
    Overwrite an existing MEICA4 output directory without waiting.

--no_para_echo
    Disable parallel echo preprocessing.
    Parallel echo preprocessing is enabled by default.
```

# Output

Assuming the command uses:

```bash
--prefix sub1_rest
```

the main outputs include:

* `sub1_rest_ctab.txt`  
  Table of component kappa, rho, variance explained, and component classifications.

* `sub1_rest_tsoc_<space>.nii.gz`  
  Optimally combined raw multi-echo time series in signal units.

* `sub1_rest_medn_<space>.nii.gz`  
  Denoised BOLD time series after basic preprocessing, T2*-weighted echo combination, and removal of definitive non-BOLD components. This is the recommended conservative denoised time-series output for many task-based activation and resting-state connectivity analyses.

* `sub1_rest_T1c_medn_<space>.nii.gz`  
  Denoised time series with T1-equilibration correction.

* `sub1_rest_hikts_<space>.nii.gz`  
  High-kappa denoised time series with T1-equilibration correction and without thermal-noise components.

* `sub1_rest_dr2s_<space>.nii.gz`  
  Denoised BOLD-like time series expressed as apparent `dR2*` units, reported in `s^-1`.

* `sub1_rest_mefc_<space>.nii.gz`  
  Accepted BOLD-like ICA coefficient maps used for ME-ICR seed-based functional connectivity analysis.

* `sub1_rest_mefl_<space>.nii.gz`  
  Full ICA coefficient image stack containing all component maps before component-selection filtering.

* `sub1_rest_mefcz_<space>.nii.gz`  
  Z-normalized spatial component maps.

* `sub1_rest_t2s_<space>.nii.gz`  
  T2*/load-related metrics produced by MEICA4.

* `sub1_rest_mmix.1D`  
  Component mixing matrix.

* `./meica.<run_label>/`  
  Directory containing preprocessing intermediate files, selected retained intermediates, and the generated processing script.

Examining accepted components in `sub1_rest_ctab.txt` together with component maps in `sub1_rest_mefl_<space>.nii.gz`, using AFNI, FSL, or another NIfTI-compatible viewer, is important for evaluating component selection quality.

Note that providing a prefix is not required, in which case the prefix of the input ME-fMRI datasets will be used to determine output file names.

## Output space suffixes

Output suffixes indicate the spatial normalization state:

* `epi`  
  Original EPI geometry.

* `nat`  
  Native anatomical space after affine anatomical-functional alignment.

* `afw`  
  Affine-warped standard space, for example affine MNI space.

* `nlw`  
  Nonlinearly warped standard space, for example using AFNI `3dQwarp`.

## Downstream analysis

MEICA4 outputs are standard NIfTI and text files suitable for use in downstream neuroimaging workflows, including AFNI, FSL, SPM, CONN, Python, R, and custom analysis pipelines. The distributed AFNI runtime is used to make MEICA4 installation and execution reproducible, but downstream analysis does not require users to remain inside an AFNI-centered workflow.

The primary denoised time-series outputs can be used for functional activation, resting-state connectivity, precision functional mapping, longitudinal analysis, and related multi-echo fMRI studies. Component maps and component tables are also exported for quality control, component review, and advanced ME-ICA/ME-ICR analyses.

# Architecture and design considerations

MEICA4 is designed around a middle ground between ease of use, reproducibility, and a controlled scientific software environment. The goal to provide a practical installation and execution model that makes multi-echo fMRI accessible to users who may not be experts in Linux system administration, AFNI installation, Python dependency management, or compiled neuroimaging software stacks.

The MEICA4 architecture supports several goals at once:

* a simple user-facing installation path;
* tight matching between MEICA4 releases and runtime releases;
* reproducible behavior across Linux and WSL2 systems;
* minimal need for administrator privileges;
* minimal global changes to the user's shell or operating system;
* minimal and controlled Python dependencies;
* predictable compiled-tool behavior for AFNI-dependent processing;
* standard NIfTI outputs that remain usable in AFNI, FSL, SPM, CONN, Python, R, and other downstream workflows.

MEICA4 is structured to make the difficult parts of ME-fMRI installation and runtime control transparent, while keeping the scientific outputs open, portable, and suitable for a wide range of downstream functional activation, connectivity, precision mapping, and longitudinal analysis workflows.

## Installer

The MEICA4 installer creates a dedicated micromamba environment for the Python layer. This gives MEICA4 a versioned Python environment without modifying the user's base Python, system Python, conda installation, or other project environments. The intent is to keep the MEICA4 release, Python dependencies, command-line behavior, and runtime assumptions tightly matched. A user should be able to install MEICA4, activate it, and run the same processing logic used by the corresponding release tag.

The distributed AFNI runtime follows the same principle for compiled neuroimaging dependencies. On Linux and WSL2, the installer downloads a runtime release matched to the MEICA4 version and installs it inside the MEICA4 source tree. `meica.py` then activates that runtime internally for its own process and for generated `_meica_*.sh` scripts. This provides tight control over the AFNI tools, template assets, OpenMP behavior, X11/Mesa/Motif-linked command-line utilities, and selected helper utilities needed by MEICA4, while avoiding broad changes to the user's shell environment.

A major architectural priority is minimizing system software requirements. Many neuroimaging tools are easy to use once installed, but installation can require system packages, administrator privileges, shell configuration, or institutional IT support. MEICA4 is intended to avoid those barriers wherever possible. The Linux/WSL2 installer should not require administrator privileges to install or use MEICA4. It should not require users to install AFNI, modify system libraries, or resolve low-level compiled dependency issues before running multi-echo fMRI analyses.

Another priority is keeping Python dependencies minimal and purposeful. MEICA4 avoids depending on numerous single-purpose Python packages. MEICA4 instead only uses widely used and tested Python scientific software, such as Numpy and Scikit-learn. This reduces installation complexity, lowers the risk of dependency conflicts, and helps avoid vulnerabilities or maintenance risks associated with useful but unmaintained packages. It also encourages performance-sensitive operations to remain in optimized compiled tools or targeted implementations within the MEICA4 source rather than relying on unoptimized Python code for heavy image-processing steps.

## Why include a distributed AFNI runtime?

Even a minimal AFNI installation often requires system-level dependencies, environment, package, and shell configuration steps. These steps commonly require administrator privileges or institutional IT support, especially on shared workstations, managed laptops, clusters, and Windows systems using WSL2. In turn, they represent barriers to entry to the ME-fMRI and ME-ICA ecosystem.

Still MEICA4 depends on AFNI for a focused set of high-quality command-line tools. The distributed runtime removes the need for users to solve the broader AFNI installation problem before they can begin working with ME-fMRI. It provides the AFNI tools, runtime libraries, templates, and selected utilities needed by MEICA4 in a controlled, reproducible package.

This approach also helps make MEICA4 outputs more accessible to downstream analysis environments such as FSL, SPM, CONN, Python-based workflows, R-based workflows, and other neuroimaging pipelines. By making the AFNI dependency transparent to the user, MEICA4 can function as a self-contained preprocessing and decomposition step whose outputs are then available for a broad range of connectivity, functional activation, precision mapping, and longitudinal studies. In that sense, the distributed runtime can make ME-fMRI analysis as accessible as conventional single-echo workflows built around FSL, SPM, or CONN, and potentially more accessible for users who want artifact-resolved multi-echo outputs without first managing AFNI-specific system dependencies.

On Linux and WSL2, the installer also downloads a matched ME-ICA AFNI runtime release and installs it inside the MEICA4 source directory:

```text
$INSTALLDIR/
  me-ica-4.0.1/
    install.sh
    meica.py
    libmeica/
    meica-afni-runtime/
      activate_meica_afni
      bin/
      lib/
      utils/
        dcm2niix
  micromamba/
  envs/
```

The full AFNI runtime is used internally by `meica.py` and by generated `_meica_*.sh` scripts. It is not globally exposed on the user's `PATH`, so it should not interfere with an existing AFNI installation. Selected user-facing utilities from the runtime, such as `dcm2niix`, are exposed after activating the MEICA4 environment.

# Additional Notes

* For more information on T2*-weighted anatomical-functional coregistration, see:  
  https://pmc.ncbi.nlm.nih.gov/articles/PMC6319659/pdf/nihms-1001520.pdf

* For more information on degrees of freedom in MEICA4-denoised data, see:  
  https://www.pnas.org/doi/abs/10.1073/pnas.1301725110

* Spatial FWHM smoothing is not recommended as a default preprocessing step. Multi-echo optimal combination already provides a tSNR benefit. For better spatial overlap of activation or connectivity patterns across subjects, use nonlinear standard-space normalization instead, for example:

```bash
meica.py ... --qwarp
```

# Citation

Please cite:

```text
MEICA4: Quantitative Precision Functional Mapping:
Stable high-dimensional fMRI decomposition and artifact-resolved denoising
of multi-echo fMRI

Prantik Kundu & R. Nathan Spreng
Whitepaper: 10.17605/OSF.IO/DN9BM
```

Additional foundational references:

```text
Kundu, P., Brenowitz, N.D., Voon, V., Worbe, Y., Vertes, P.E., Inati, S.J.,
Saad, Z.S., Bandettini, P.A. & Bullmore, E.T. Integrated strategy for improving
functional connectivity mapping using multiecho fMRI. PNAS (2013).

Kundu, P., Inati, S.J., Evans, J.W., Luh, W.M. & Bandettini, P.A.
Differentiating BOLD and non-BOLD signals in fMRI time series using multi-echo EPI.
NeuroImage (2011).
```
