Introduction
============

``meica`` preprocesses multi-echo datasets and applies multi-echo ICA based
on spatially concatenated echoes. It does so in the following steps:

#. Calculates motion parameters based on images with highest contrast (usually the first echo)
#. Applies motion correction and T2*-weighted co-registration parameters
#. Applies standard EPI preprocessing (slice-time correction, etc.)
#. Computes PCA and ICA in conjunction with TE-dependence analysis

Derivatives
-----------

* ``*_medn.nii.gz``
    'Denoised' BOLD time series after: basic preprocessing,
    T2* weighted averaging of echoes (i.e. 'optimal combination'),
    ICA denoising.
    Use this dataset for task analysis and resting state time series correlation analysis.
* ``*_tsoc.nii.gz``
    'Raw' BOLD time series dataset after: basic preprocessing
    and T2* weighted averaging of echoes (i.e. 'optimal combination').
    'Standard' denoising or task analyses can be assessed on this dataset
    (e.g. motion regression, physio correction, scrubbing, etc.)
    for comparison to ME-ICA denoising.
* ``*_mefc.nii.gz``
    Component maps (in units of \delta S) of accepted BOLD ICA components.
    Use this dataset for ME-ICR seed-based connectivity analysis.
* ``*_mefl.nii.gz``
    Component maps (in units of \delta S) of ALL ICA components.
* ``*_ctab.nii.gz``
    Table of component Kappa, Rho, and variance explained values, plus listing of component classifications.
