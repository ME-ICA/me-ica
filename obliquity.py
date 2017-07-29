#!/usr/bin/env python

import nibabel as nib
import numpy as np

def check_obliquity(dset):
    """
    Determines if `dset` is oblique

    Parameters
    ----------
    dset : str
        path to file

    Returns
    -------
    bool : whether `dset` is oblique (True)
    """

    aff = nib.load(dset).affine

    dxtmp = np.sqrt((aff[:3,0]**2).sum())
    xmax = np.abs(aff[:3,0]).max() / dxtmp

    dytmp = np.sqrt((aff[:3,1]**2).sum())
    ymax = np.abs(aff[:3,1]).max() / dytmp

    dztmp = np.sqrt((aff[:3,2]**2).sum())
    zmax = np.abs(aff[:3,2]).max() / dztmp

    fig_merit = np.min([xmax, ymax, zmax])
    ang_merit = (np.arccos(fig_merit) * 180) / np.pi

    return ang_merit != 0.0
