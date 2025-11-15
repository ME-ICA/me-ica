from typing import Union

import nibabel as nib
import numpy as np

import scipy.stats as stats

from .selection import scoreatpercentile

nhs = np.hstack
nvs = np.vstack


def niwrite(data, affine, name, header):
    print(
        " + Writing file: %s ...." % name,
    )

    head = header
    thishead = head.copy()
    thishead.set_data_shape(list(data.shape))

    outni = nib.Nifti1Image(data, affine, header=thishead)  # type: ignore
    outni.to_filename(name)
    print("done.")


def makemask(cdat):
    nx, ny, nz, Ne, nt = cdat.shape

    mask = np.ones((nx, ny, nz), dtype=bool)

    for i in range(Ne):
        tmpmask = (cdat[:, :, :, i, :] != 0).prod(axis=-1, dtype=bool)
        mask = mask & tmpmask

    return mask


def fmask(data, mask):
    """
    fmask(data,mask)

    Input:
    data shape is (nx,ny,nz,...)
    mask shape is (nx,ny,nz)

    Output:
    out shape is (Nm,...)
    """

    s = data.shape

    N = s[0] * s[1] * s[2]
    news = []
    news.append(N)

    if len(s) > 3:
        news.extend(s[3:])

    tmp1 = np.reshape(data, news)
    fdata = tmp1.compress((mask > 0).ravel(), axis=0)

    return fdata.squeeze()


def cat2echos(data, Ne):
    """
    cat2echos(data,Ne)

    Input:
    data shape is (nx,ny,Ne*nz,nt)
    """
    nx, ny = data.shape[0:2]
    nz = data.shape[2] // Ne
    if len(data.shape) > 3:
        nt = data.shape[3]
    else:
        nt = 1
    return np.reshape(data, (nx, ny, nz, Ne, nt), order="F")


def uncat2echos(data, Ne):
    """
    uncat2echos(data,Ne)

    Input:
    data shape is (nx,ny,Ne,nz,nt)
    """
    nx, ny = data.shape[0:2]
    nz = data.shape[2] * Ne
    if len(data.shape) > 4:
        nt = data.shape[4]
    else:
        nt = 1
    return np.reshape(data, (nx, ny, nz, nt), order="F")


def unmask(data, mask):
    """
    unmask (data,mask)

    Input:

    data has shape (Nm,nt)
    mask has shape (nx,ny,nz)

    """
    M = (mask != 0).ravel()
    Nm = M.sum()

    nx, ny, nz = mask.shape

    if len(data.shape) > 1:
        nt = data.shape[1]
    else:
        nt = 1

    out = np.zeros((nx * ny * nz, nt), dtype=data.dtype)
    out[M, :] = np.reshape(data, (Nm, nt))

    return np.squeeze(np.reshape(out, (nx, ny, nz, nt)))


def eimask(dd, ees=None):
    if ees is None:
        ees = list(range(dd.shape[1]))
    imask = np.zeros([dd.shape[0], len(ees)])
    for ee in ees:
        print(ee)
        lthr = 0.001 * scoreatpercentile(dd[:, ee, :].flatten(), 98)
        hthr = 5 * scoreatpercentile(dd[:, ee, :].flatten(), 98)
        print(lthr, hthr)
        imask[dd[:, ee, :].mean(1) > lthr, ee] = 1
        imask[dd[:, ee, :].mean(1) > hthr, ee] = 0
    return imask


def makeadmask(cdat, min=True, getsum=False):
    nx, ny, nz, Ne, nt = cdat.shape

    mask = np.ones((nx, ny, nz), dtype=bool)

    if min:
        mask = cdat[:, :, :, :, :].prod(axis=-1).prod(-1) != 0
        return mask
    else:
        # Make a map of longest echo that a voxel can be sampled with,
        # with minimum value of map as X value of voxel that has median
        # value in the 1st echo. N.b. larger factor leads to bias to lower TEs
        emeans = cdat[:, :, :, :, :].mean(-1)
        medv = emeans[:, :, :, 0] == stats.scoreatpercentile(
            emeans[:, :, :, 0][emeans[:, :, :, 0] != 0],
            33,
            interpolation_method="higher",
        )
        lthrs = np.squeeze(
            np.array([emeans[:, :, :, ee][medv] / 3 for ee in range(Ne)])
        )
        if len(lthrs.shape) == 1:
            lthrs = np.atleast_2d(lthrs).T
        lthrs = lthrs[:, lthrs.sum(0).argmax()]
        lthrs[lthrs < 0] = 0
        mthr = np.ones([nx, ny, nz, Ne])
        for ee in range(Ne):
            mthr[:, :, :, ee] *= lthrs[ee]
        mthr = np.abs(emeans[:, :, :, :]) > mthr
        masksum = np.array(mthr, dtype=int).sum(-1)
        mask = masksum != 0
        if getsum:
            return mask, masksum
        else:
            return mask


def is_tedlike(tedarr: np.ndarray) -> bool:
    if len(tedarr.shape) == 4 and tedarr.shape[3] >= 9:  # TODO: Change 9 to Constant
        return True
    return False


def place3d(
    src_box,
    src_anchor,
    tgt_box,
    tgt_anchor,
    *,
    sourceval: Union[int, float, None] = None,
    inplace=False,
    add=False,
):
    assert np.min(src_anchor) >= 0
    assert np.min(tgt_anchor) >= 0

    src_anchor = np.round(src_anchor)
    tgt_anchor = np.round(tgt_anchor)

    tgt_box_shape = tgt_box.shape
    if not inplace and not add:
        tgt_box = np.zeros(tgt_box_shape)

    src_shape = src_box.shape

    start_source = -src_anchor  # think of src_anchor like cmass
    end_source = start_source + src_shape

    start_target = tgt_anchor - src_anchor
    # tgt_anchor + src_shape - src_anchor, switching order for type sanity
    end_target = -src_anchor + tgt_anchor + src_shape

    adj_start_target = np.max([[0, 0, 0], start_target], axis=0)
    adj_end_target = np.min([tgt_box_shape, end_target], axis=0)

    start_adj = adj_start_target - start_target
    end_adj = adj_end_target - end_target

    start_source_adj = start_source + start_adj
    end_source_adj = end_source + end_adj

    start_source_abs = start_source_adj + src_anchor
    end_source_abs = end_source_adj + src_anchor

    _sp = start_source_abs.astype(int)
    _ep = end_source_abs.astype(int)
    _sr = adj_start_target.astype(int)
    _er = adj_end_target.astype(int)

    if sourceval is not None:
        _src_box = np.array(src_box != 0)
        _src_box[src_box != 0] = sourceval
        src_box = _src_box

    toplace = src_box[_sp[0] : _ep[0], _sp[1] : _ep[1], _sp[2] : _ep[2]]  # noqa # type: ignore

    if add:
        tgt_box[_sr[0] : _er[0], _sr[1] : _er[1], _sr[2] : _er[2]] = (  # noqa
            tgt_box[_sr[0] : _er[0], _sr[1] : _er[1], _sr[2] : _er[2]] + toplace  # noqa
        )
    else:
        tgt_box[_sr[0] : _er[0], _sr[1] : _er[1], _sr[2] : _er[2]] = toplace  # noqa

    return tgt_box
