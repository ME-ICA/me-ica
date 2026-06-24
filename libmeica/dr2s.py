import nibabel as nib
import numpy as np

from libmeica.utils.selection import z_
from libmeica.utils.volume import niwrite, ts_sigma_rescale, unmask

__version__ = "4.0.1"


def signal_to_dr2s(*, assets, glsig=None, lambda_r2star=0):
    """
    Convert dataset from signal units to t2s units
    For now we're expecting datasets to have been written
    we're reloading them and doing this calculation.
    There's an in-memory way to do this but this
    implementation is good for proof-of-concept
    and testing.
    """
    tes = np.array(assets.tes) * 0.001  # in seconds
    hik_fn = None
    if glsig is None:
        ts_fns = [f"hik_ts_e{_i}.nii" for _i in range(1, len(tes) + 1)]
        hik_fn = "hik_ts_OC.nii"
    else:
        ts_fns = [f"hik_ts_e{_i}_T1c.nii" for _i in range(1, len(tes) + 1)]
        hik_fn = "hik_ts_OC_T1c.nii"  # added 9/19

    t2sG = assets.t2sG
    t2sGm = t2sG != 0
    t2sm = t2sG[t2sGm]
    mum = assets.mu[t2sGm]

    # Read in split time series in delta-S
    ts_vols = [nib.load(_ts_fn) for _ts_fn in ts_fns]  # type: ignore
    ts_dats = [_ts_vol.get_fdata() for _ts_vol in ts_vols]  # type: ignore
    ts_sigmas = [_ts_dat[t2sGm].std(1) for _ts_dat in ts_dats]
    B = np.array(ts_sigmas)  # * np.sqrt(nT)

    # t2s keep in volume
    NmD = t2sm.flatten().shape[0]
    X2 = (-1 * np.tile(np.atleast_2d(tes).T, (1, NmD)).T * mum).T

    # Fit for dR2*
    r2coef = -1 * (B * X2).sum(axis=0) / (X2**2).sum(axis=0)  # dR2*
    r2coef_z = z_(r2coef)
    r2coef[r2coef_z > 6] = 0
    r2coef[t2sG[t2sGm] > assets.T2brain] = 0

    # Write out amplitudes as delta-R2*
    niwrite(unmask(r2coef, t2sGm), assets.aff, "ampdr2.nii", header=assets.head)

    ts_sigma_rescale(
        hik_fn, scale_fac=unmask(r2coef, assets.mask), out_suffix=".dr2s.nii"
    )
