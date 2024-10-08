import numpy as np

from .volume import fmask, niwrite, unmask
from .selection import scoreatpercentile


def t2smap(catd, mask, tes):
    """
    t2smap(catd,mask,tes)

    Input:

    catd  has shape (nx,ny,nz,Ne,nt)
    mask  has shape (nx,ny,nz)
    tes   is a 1d numpy array
    """
    nx, ny, nz, Ne, nt = catd.shape
    N = nx * ny * nz

    echodata = fmask(catd, mask)
    Nm = echodata.shape[0]

    # Do Log Linear fit
    B = np.reshape(np.abs(echodata[:, :Ne]) + 1, (Nm, (Ne) * nt)).transpose()
    B = np.log(B)
    x = np.array([np.ones(Ne), -tes])
    X = np.tile(x, (1, nt))
    X = np.sort(X)[:, ::-1].transpose()

    beta, res, rank, sing = np.linalg.lstsq(X, B)
    t2s = 1 / beta[1, :].transpose()
    s0 = np.exp(beta[0, :]).transpose()

    # Goodness of fit
    alpha = (np.abs(B) ** 2).sum(axis=0)
    t2s_fit = (alpha - res) / (2 * res)

    out = (
        np.squeeze(unmask(t2s, mask)),
        np.squeeze(unmask(s0, mask)),
        unmask(t2s_fit, mask),
    )

    return out


def t2sadmap(catd, mask, tes, masksum, assets=None):
    """
    t2smap(catd,mask,tes)

    Input:

    catd  has shape (nx,ny,nz,Ne,nt)
    mask  has shape (nx,ny,nz)
    tes   is a 1d numpy array
    """
    nx, ny, nz, Ne, nt = catd.shape
    N = nx * ny * nz

    echodata = fmask(catd, mask)
    Nm = echodata.shape[0]

    t2ss = np.zeros([nx, ny, nz, Ne - 1])
    s0s = np.zeros([nx, ny, nz, Ne - 1])

    for ne in range(1, Ne + 1):
        # Do Log Linear fit
        B = np.reshape(np.abs(echodata[:, :ne]) + 1, (Nm, (ne) * nt)).transpose()
        B = np.log(B)
        x = np.array([np.ones(ne), -tes[:ne]])
        X = np.tile(x, (1, nt))
        X = np.sort(X)[:, ::-1].transpose()

        beta, res, rank, sing = np.linalg.lstsq(X, B, rcond=None)
        t2s = 1 / beta[1, :].transpose()
        s0 = np.exp(beta[0, :]).transpose()

        t2s[np.isinf(t2s)] = 500.0
        s0[np.isnan(s0)] = 0.0

        t2ss[:, :, :, ne - 2] = np.squeeze(unmask(t2s, mask))
        s0s[:, :, :, ne - 2] = np.squeeze(unmask(s0, mask))

    # Limited T2* and S0 maps
    fl = np.zeros([nx, ny, nz, len(tes) - 2 + 1])
    for ne in range(Ne - 1):
        fl_ = np.squeeze(fl[:, :, :, ne])
        fl_[masksum == ne + 2] = True
        fl[:, :, :, ne] = fl_
    fl = np.array(fl, dtype=bool)
    t2s = np.squeeze(unmask(t2ss[fl], masksum > 1))
    s0 = np.squeeze(unmask(s0s[fl], masksum > 1))

    # Full T2* maps with S0 estimation errors
    t2sG = t2s.copy()
    s0G = s0.copy()
    t2sG[masksum == 1] = t2ss[masksum == 1, 0]
    s0G[masksum == 1] = s0s[masksum == 1, 0]

    if assets is not None:
        # Condition values
        cap_t2s = scoreatpercentile(t2s.flatten(), 99.5)
        t2s[t2s > cap_t2s * 10] = cap_t2s
        niwrite(t2s, assets.aff, "t2sv.nii", header=assets.head)
        niwrite(s0, assets.aff, "s0v.nii", header=assets.head)
        niwrite(t2ss, assets.aff, "t2ss.nii", header=assets.head)
        niwrite(s0s, assets.aff, "s0vs.nii", header=assets.head)
        niwrite(t2sG, assets.aff, "t2svG.nii", header=assets.head)
        niwrite(s0G, assets.aff, "s0vaf.nii", header=assets.head)

    return t2s, s0, t2ss, s0s, t2sG, s0G
