"""
spatialfilters.py: A set of functions to process volumes for spatial
   maniputations

Author: Prantik Kundu, PhD
"""

import datetime
import os
import subprocess
from typing import Union
from uuid import uuid4

import nibabel as nib
import numpy as np

# import scipy.spatial.distance as distance
from scipy import ndimage

from .volume import fmask, unmask

# _a = np.array


def niwrite(data, affine, name, header=None, quiet=False):
    """
    Temporary utility function. Main function will be moved to
    IO module with lib_tedana refactoring
    """
    # TO DO: push niwrite out to an io module
    if not quiet:
        print(" + Writing file: %s ...." % name)

    # thishead = headerV
    # if thishead == None:
    #     thishead = head.copy()
    #     thishead.set_data_shape(list(data.shape))
    outni = nib.Nifti1Image(data, affine, header=header)  # type: ignore
    outni.to_filename(name)
    if not quiet:
        print("done.")


def localfwhm_afni(
    data,
    maskvol,
    *,
    header,
    affine,
    nhtype="SPHERE",
    nhsize=-1.42,  # Indicating 19 pixel hood for FWHM
    thr=0,
    domask=True,
    inpath="",
    inbrik=0,
    maskpath="",
    prefix="",
    outdir="/tmp/",
    rmoutfile=True,
):
    rminfile = rmmaskfile = False
    uid = uuid4()  # datetime.datetime.now().isoformat()

    # Write the input file
    if inpath == "":
        if data.ndim == 4 or data.ndim == 2:
            data = data[..., inbrik]
        else:
            data = data.copy()
        assert thr >= 0
        data[np.abs(data) < thr] = 0
        if np.squeeze(data).ndim == 1:
            _data = unmask(data, maskvol)
        else:
            _data = data
        _inprefix = f"__lsin_{uid}"
        _inpath = f"{outdir}/{_inprefix}.nii.gz"
        niwrite(_data, affine, _inpath, header, quiet=True)
        rminfile = True
    else:
        _inpath = f"{inpath}"

    # Write the mask file
    if maskpath == "":
        _maskpath = f"{outdir}/_{uid}_localfwhm_mask.nii.gz"
        niwrite(maskvol, affine, _maskpath, header, quiet=True)
        rmmaskfile = True
    else:
        _maskpath = maskpath

    # Define the prefix
    if prefix == "":
        _outpath = f"{_inpath}_localfwhm_{uid}.nii.gz"
    else:
        _outpath = f"{outdir}/{prefix}.nii.gz"

    # Run the command
    if data is not None:
        _inbrik = 0
    else:
        _inbrik = inbrik
    runcmd = f"3dLocalstat -overwrite -nbhd '{nhtype}({nhsize})' -stat FWHM -mask {_maskpath} -prefix {_outpath} {_inpath}[{_inbrik}]"
    run_result = subprocess.run(
        runcmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
    )

    # Read in
    analyzed = nib.load(_outpath).get_fdata()  # type: ignore
    analyzed[~maskvol] = 0
    if domask:
        analyzed = fmask(analyzed, maskvol)

    # Allow deletion of files before return
    analyzed = np.asarray(np.squeeze(analyzed))

    # Clean up
    if rminfile:
        os.system(f"rm {_inpath}")
    if rmmaskfile:
        os.system(f"rm {_maskpath}")
    if rmoutfile:
        os.system(f"rm {_outpath}")

    return analyzed


def threshold(voltothr, thr):
    thresholded = np.zeros(voltothr.shape)
    thresholded[voltothr > thr] = voltothr[voltothr > thr]
    return thresholded


def thr_z(arr, z=1.0, f=1.0):
    _mn = np.mean(arr)
    _std = np.std(arr)
    _thr = _mn * f + z * _std
    mask = np.abs(arr) > _thr
    return mask


def spatclust(arr_in, thr, csize, *, mask=None, getsizes=False):
    """
    Uses scipy.ndimage.label for clustering

    TODO: Type annotations!

    Input:
        data: 3-D numpy array of image
        thr: float
        csize: int
    Output:
        label_mask: 3-D numpy array of clustered and
            labelled data
    """

    # Make sure we're dealing with a volume with positive values
    if np.sum(arr_in < 0) != 0:
        # print("spatclust() does not support negative values. Zeroing negatives")
        arr_in = arr_in.copy()
        arr_in[arr_in < 0] = 0

    # Apply a volume mask if provided
    if mask is not None:
        arr = unmask(arr_in, mask)
    else:
        arr = arr_in

    assert arr.ndim == 3

    arr_shape = arr.shape
    arr = np.squeeze(arr)
    squeeze_ndim = arr.ndim
    struc = ndimage.generate_binary_structure(squeeze_ndim, 1)

    # Initialize output
    clustered_out = np.zeros(arr.shape, int)

    # Make a binary mask
    binmask = arr > thr

    # Cluster the data with labels
    labeled, _ = ndimage.label(binmask, struc)  # type: ignore

    # Fast tally of clusters
    unique, counts = np.unique(labeled, return_counts=True)
    clust_sizes = dict(zip(unique, counts))
    clust_sizes = {k: v for k, v in clust_sizes.items() if v >= csize}

    # Write outputs
    for clust_i, clust_l in enumerate(clust_sizes.keys()):
        if np.all(binmask[labeled == clust_l] == 1):
            clustered_out[labeled == clust_l] = clust_i

    # Remask
    if mask is not None:
        arr_return = fmask(clustered_out, mask)
    else:
        arr_return = clustered_out

    # Recast singlteton dim
    if squeeze_ndim != 3:
        arr_return = np.reshape(arr_return, arr_shape)

    if getsizes:
        return arr_return, clust_sizes
    else:
        return arr_return


def spatclust_bi(arr_in, thr, csize, mask=None, both=False):
    """
    Conducts spatial clustering on volumes with positives and negatives,
    clustered separately or together (abs)

    Input:
        data: 3-D numpy array of image
        thr: float
        csize: int
        both: bool
    Output:
        clust_tot: Total 3-D numpy array of clustered and
            labelled data for positive and negative values
    """

    if both:
        clust_tot = spatclust(np.abs(arr_in), thr, csize, mask=mask)

    else:
        clust_pos = spatclust(arr_in, thr, csize, mask=mask)
        clust_neg = spatclust(-arr_in, thr, csize, mask=mask)
        neg_max = np.max(clust_neg)  # type: ignore
        clust_tot = clust_neg + clust_pos + neg_max * (clust_pos != 0)

    return clust_tot


def gradmask(img, norm=True):
    # Get a treshold for psc based on that at a high z-value?
    #   nope z-vlaues correspond to p-s that are too low
    # Calculate threshold params
    nbrvox = (img != 0).sum()
    clust_size = int(0.0005 * nbrvox)
    psc_val_mask = thr_z(img, z=8.0)
    img_grad = np.gradient(img)
    img_grad_norm = np.linalg.norm(img_grad, axis=0)
    img_grad_mask = thr_z(img_grad_norm, z=0.5, f=1)
    grad_mask = psc_val_mask ^ img_grad_mask
    grad_cl = spatclust(grad_mask, csize=clust_size, thr=0)
    _struc = ndimage.generate_binary_structure(3, 3)
    grad_cl = ndimage.binary_closing(grad_cl, structure=_struc)
    if norm:
        grad_cl = grad_cl * img
    else:
        grad_cl = grad_cl[..., np.newaxis] * np.array(img_grad).transpose(1, 2, 3, 0)
    return grad_cl


def gradmask_sep(img, z=1.0, clfac=0.0005, closingstep=3):
    """
    Compute a clustered mask of high gradients separately
    for each direction
    """
    # Set clustering parameters
    nbrvox = (img != 0).sum()
    clust_size = int(clfac * nbrvox)
    closing_struc = ndimage.generate_binary_structure(3, closingstep)

    # Thresholded grad directions
    img_grad = np.gradient(img)
    thr_grad = []
    for _dirimg in img_grad:
        thr_dir = thr_z(_dirimg, z) * _dirimg
        grad_cl = spatclust(thr_dir, csize=clust_size, thr=0)
        grad_cl = ndimage.binary_closing(grad_cl, structure=closing_struc)
        thr_grad.append(grad_cl * thr_dir)

    return np.array(thr_grad).transpose(1, 2, 3, 0)


def rankvec(vals):
    asort = np.argsort(vals)
    ranks = np.zeros(vals.shape[0])
    ranks[asort] = np.arange(vals.shape[0]) + 1
    return ranks
