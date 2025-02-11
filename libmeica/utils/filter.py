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


def niwrite(data, affine, name, header=None):
    """
    Temporary utility function. Main function will be moved to
    IO module with lib_tedana refactoring
    """
    # TO DO: push niwrite out to an io module
    print(" + Writing file: %s ...." % name)

    # thishead = headerV
    # if thishead == None:
    #     thishead = head.copy()
    #     thishead.set_data_shape(list(data.shape))
    outni = nib.Nifti1Image(data, affine, header=header)  # type: ignore
    outni.to_filename(name)
    print("done.")


def spatclust_afni(
    data, mask, csize, thr, header, aff, infile=None, dindex=0, tindex=0, prefix=""
):
    """
    Deprecated
    """

    if infile is None:
        data = data.copy()
        data[data < thr] = 0
        niwrite(unmask(data, mask), aff, "__clin.nii.gz", header)
        infile = "__clin.nii.gz"
    addopts = ""
    if data is not None and len(np.squeeze(data).shape) > 1 and dindex + tindex == 0:
        addopts = "-doall"
    else:
        addopts = "-1dindex %s -1tindex %s" % (str(dindex), str(tindex))
    os.system(
        "3dmerge -overwrite %s -dxyz=1  -1clust 1 %i -1thresh %.02f -prefix %s__clout.nii.gz %s"
        % (addopts, int(csize), float(thr), prefix, infile)
    )
    clustered = (
        fmask(nib.load("%s__clout.nii.gz" % prefix).get_fdata(), mask) != 0  # type: ignore
    )  # changed here too
    os.system("rm %s__clout.nii.gz" % prefix)
    return clustered


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
        niwrite(_data, affine, _inpath, header)
        rminfile = False
    else:
        _inpath = f"{inpath}"

    # Write the mask file
    if maskpath == "":
        _maskpath = f"{outdir}/_{uid}_localfwhm_mask.nii.gz"
        niwrite(maskvol, affine, _maskpath, header)
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
    run_result = subprocess.run(runcmd, shell=True)

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
        print("spatclust() does not support negative values. Zeroing negatives")
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


# def fmask(data, mask):
#     """
#     fmask(data,mask)

#     Input:
#     data shape is (nx,ny,nz,...)
#     mask shape is (nx,ny,nz)

#     Output:
#     out shape is (Nm,...)
#     """

#     s = data.shape
#     # sm = mask.shape

#     N = s[0] * s[1] * s[2]
#     news = []
#     news.append(N)

#     if len(s) > 3:
#         news.extend(s[3:])

#     tmp1 = np.reshape(data, news)
#     fdata = tmp1.compress((mask > 0).ravel(), axis=0)

#     return fdata.squeeze()


# def unmask(data, mask):
#     """
#     unmask (data,mask)

#     Input:

#     data has shape (Nm,nt)
#     mask has shape (nx,ny,nz)

#     """
#     M = (mask != 0).ravel()
#     Nm = M.sum()

#     nx, ny, nz = mask.shape

#     if len(data.shape) > 1:
#         nt = data.shape[1]
#     else:
#         nt = 1

#     out = np.zeros((nx * ny * nz, nt), dtype=data.dtype)
#     out[M, :] = np.reshape(data, (Nm, nt))

#     return np.squeeze(np.reshape(out, (nx, ny, nz, nt)))


# def com3d(arr):
#     """
#     Find the center of mass of an array. Can return an empty
#     COM given an irregular object. Use com3d_enclosed to prevent
#     against that.
#     """
#     assert arr.ndim == 3

#     # Calculate the total mass (sum of all values in the 3D array)
#     total_mass = arr.sum()
#     assert total_mass != 0

#     x_sum = (np.arange(arr.shape[0])[:, np.newaxis, np.newaxis] * arr).sum()
#     y_sum = (np.arange(arr.shape[1])[np.newaxis, :, np.newaxis] * arr).sum()
#     z_sum = (np.arange(arr.shape[2])[np.newaxis, np.newaxis, :] * arr).sum()

#     # Calculate the center of mass
#     com_x = x_sum / total_mass
#     com_y = y_sum / total_mass
#     com_z = z_sum / total_mass

#     return np.array((com_x, com_y, com_z))


# def com3d_enclosed(arr):
#     """
#     Returns center of mass or closest enclosed point of the
#     target object in arr
#     """
#     arr = arr.astype(int)
#     com_cand_coor = tuple(com3d(arr).astype(int))
#     com_cand = com3d(arr).astype(int)
#     if arr.item(com_cand_coor) == 0:
#         _idx = np.indices(arr.shape)
#         nz_idx = _idx[..., arr != 0]
#         # [..., arr != 0]
#         dists = np.linalg.norm(nz_idx - com_cand[:, np.newaxis], axis=0)
#         com_enc = nz_idx[:, np.argmin(dists)]
#         com_enc_coor = tuple(com_enc.astype(int))
#         assert arr.item(com_enc_coor) != 0
#         return com_enc_coor
#     else:
#         return com_cand


# def autobox(patch: np.ndarray):
#     idx = np.indices(patch.shape)[:, patch != 0]
#     idx = idx.reshape(3, -1)
#     rng = list(zip(np.min(idx, 1), np.max(idx, 1) + 1))  # We're counting the max
#     # index here, so when using
#     # for range we add 1
#     take0 = np.take(patch, np.arange(*rng[0]), axis=0)
#     take1 = np.take(take0, np.arange(*rng[1]), axis=1)
#     take2 = np.take(take1, np.arange(*rng[2]), axis=2)
#     return take2


# def autobox_mask(patch: np.ndarray):
#     idx = np.indices(patch.shape)[:, patch != 0]
#     idx = idx.reshape(patch.ndim, -1)
#     rng = list(zip(np.min(idx, 1), np.max(idx, 1) + 1))  # We're counting the max
#     # index here, so when using
#     # for range we add 1
#     patch_map = np.zeros(patch.shape)
#     # Generate all combinations of indices using np.meshgrid
#     indices = np.meshgrid(*[np.arange(start, end) for start, end in rng], indexing="ij")
#     # Flatten the indices and update the corresponding elements to 1
#     patch_map[tuple(idx.flatten() for idx in indices)] = 1
#     return patch_map


# def comgrid3d(arr1: np.ndarray, arr2: np.ndarray, do_autobox=True):
#     """
#     Place 2 3D arrays or patches on a common grid with
#     each's COM at origin
#     """
#     # Make the common grid
#     if do_autobox:
#         ab_arr1 = autobox(arr1)
#         ab_arr2 = autobox(arr2)
#     else:
#         ab_arr1 = arr1.copy()
#         ab_arr2 = arr2.copy()
#     ab_shapes = list(zip(ab_arr1.shape, ab_arr2.shape))
#     tgt_shape = np.max(ab_shapes, 1)

#     # Determine center points of patches
#     tgt_center = tgt_shape // 2
#     com1 = np.array(com3d(ab_arr1), dtype=int)
#     com2 = np.array(com3d(ab_arr2), dtype=int)

#     # Find start and end points relative to com and autobox
#     #  dimensions
#     tgt1_start = tgt_center - com1
#     tgt2_start = tgt_center - com2
#     tgt1_end = tgt1_start + ab_arr1.shape
#     tgt2_end = tgt2_start + ab_arr2.shape

#     # Expand shapes due to out of window shifts
#     #  in center of mass
#     tgt_shape_old = tgt_shape
#     start_expand = np.min([tgt1_start, tgt2_start, [0, 0, 0]], axis=0)
#     end_expand = np.max([tgt1_end - tgt_shape, tgt2_end - tgt_shape, [0, 0, 0]], axis=0)
#     tgt_shape = tgt_shape + np.abs(start_expand) + np.abs(end_expand)
#     tgt1_start = tgt1_start - start_expand
#     tgt2_start = tgt2_start - start_expand
#     tgt1_end = tgt1_start + np.array(ab_arr1.shape)
#     tgt2_end = tgt2_start + np.array(ab_arr2.shape)

#     tgt1 = np.zeros(tgt_shape)
#     tgt2 = np.zeros(tgt_shape)

#     # Regrid input arrays into common grid
#     tgt1_rng = list(zip(tgt1_start, tgt1_end))
#     tgt2_rng = list(zip(tgt2_start, tgt2_end))

#     tgt1[
#         tgt1_rng[0][0] : tgt1_rng[0][1],  # noqa
#         tgt1_rng[1][0] : tgt1_rng[1][1],  # noqa
#         tgt1_rng[2][0] : tgt1_rng[2][1],  # noqa
#     ] = ab_arr1
#     tgt2[
#         tgt2_rng[0][0] : tgt2_rng[0][1],  # noqa
#         tgt2_rng[1][0] : tgt2_rng[1][1],  # noqa
#         tgt2_rng[2][0] : tgt2_rng[2][1],  # noqa
#     ] = ab_arr2

#     return tgt1, tgt2


# def dice3d(arr1: np.ndarray, arr2: np.ndarray, do_autobox=True) -> float:
#     """Compute Dice metric of 2 arbitrarily shaped blobs"""
#     arr1_c, arr2_c = comgrid3d(arr1, arr2, do_autobox=do_autobox)
#     # default return on dice is dissimilar
#     #  we do 1- to get similarity
#     return 1 - distance.dice(arr1_c.flatten(), arr2_c.flatten())
