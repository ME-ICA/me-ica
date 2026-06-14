import numpy as np

from libmeica.utils.memory import create_memmap

from .utils.volume import niwrite, unmask


def gscontrol_raw(*, OCcatd, assets, dtrank=4):
    """
    This function uses the spatial global signal estimation approach to modify catd (global variable) to
    removal global signal out of individual echo time series datasets. The spatial global signal is estimated
    from the optimally combined data after detrending with a Legendre polynomial basis of order=0 and degree=dtrank.
    """
    print("++ Applying amplitude-based T1 equilibration correction")

    ne = assets.Ne
    head = assets.head
    aff = assets.aff
    catd = assets.catd
    nx, ny, nz, nt = OCcatd.shape
    out_catd = np.zeros([nx, ny, nz, ne, nt])

    # Write out original OC signal
    niwrite(OCcatd, aff, "tsoc_orig.nii", head)

    # Legendre polynomial basis for denoising
    from scipy.special import lpmv

    Lmix = np.array(
        [lpmv(0, vv, np.linspace(-1, 1, OCcatd.shape[-1])) for vv in range(dtrank)]
    ).T

    # Compute mean, std, mask local to this function - inefficient, but makes this function a bit more modular
    Gmu = OCcatd.mean(-1)
    # Gstd = OCcatd.std(-1)
    Gmask = Gmu != 0

    # Find spatial global signal
    dat = OCcatd[Gmask] - Gmu[Gmask][:, np.newaxis]
    sol = np.linalg.lstsq(Lmix, dat.T, rcond=None)  # Legendre basis for detrending
    detr = dat - np.dot(sol[0].T, Lmix.T)[0]
    sphis = (detr).min(1)
    sphis -= sphis.mean()
    T1mi = unmask(sphis, Gmask)
    assets.t1mi = T1mi
    niwrite(T1mi, aff, "T1gs.nii", head)

    # Find time course of the spatial global signal, make basis with the Legendre basis
    glsig = np.linalg.lstsq(np.atleast_2d(sphis).T, dat, rcond=None)[0]
    glsig = (glsig - glsig.mean()) / glsig.std()
    np.savetxt("T1_dephase.1D", glsig)
    glbase = np.hstack([Lmix, glsig.T])

    # Project global signal out of optimally combined data
    sol = np.linalg.lstsq(np.atleast_2d(glbase), dat.T, rcond=None)
    tsoc_nogs = (
        dat
        - np.dot(np.atleast_2d(sol[0][dtrank]).T, np.atleast_2d(glbase.T[dtrank]))
        + Gmu[Gmask][:, np.newaxis]
    )

    # Overwrite OCcatd in assets, to be used in rest of workflow
    OCcatd_mir = unmask(tsoc_nogs, Gmask)
    assets.OCcatd = create_memmap("_OCcatd_mir", assets.tsshape)
    assets.OCcatd[:] = OCcatd_mir
    niwrite(OCcatd_mir, aff, "tsoc_nogs.nii", head)

    # breakpoint()

    # Project glbase out of each echo
    for ii in range(ne):
        dat = catd[:, :, :, ii, :][Gmask]
        sol = np.linalg.lstsq(np.atleast_2d(glbase), dat.T, rcond=None)
        e_nogs = dat - np.dot(
            np.atleast_2d(sol[0][dtrank]).T, np.atleast_2d(glbase.T[dtrank])
        )
        out_catd[:, :, :, ii, :] = unmask(e_nogs, Gmask)

    return out_catd


def gscontrol_mmix(OCcatd, mmix, mask, acc, assets):
    head = assets.head
    aff = assets.aff

    Gmu = OCcatd.mean(-1)
    Gstd = OCcatd.std(-1)
    Gmask = Gmu != 0

    """
    Compute temporal regression
    """
    dat = (OCcatd[Gmask] - Gmu[Gmask][:, np.newaxis]) / Gstd[mask][:, np.newaxis]
    solG = np.linalg.lstsq(mmix, dat.T, rcond=None)
    resid = dat - np.dot(solG[0].T, mmix.T)

    """
    Build BOLD time series without amplitudes, and save T1-like effect
    """
    bold_ts = np.dot(solG[0].T[:, acc], mmix[:, acc].T)
    sphis = bold_ts.min(-1)
    sphis -= sphis.mean()
    print(sphis.shape)
    niwrite(unmask(sphis, mask), aff, "sphis_hik.nii", header=head)

    """
    Find the global signal based on the T1-like effect
    """
    sol = np.linalg.lstsq(np.atleast_2d(sphis).T, dat, rcond=None)
    glsig = sol[0]
    np.savetxt("T2_dephase.1D", glsig)

    """
    T1 correct time series by regression
    """
    bold_noT1gs = bold_ts - np.dot(
        np.linalg.lstsq(glsig.T, bold_ts.T, rcond=None)[0].T, glsig
    )
    niwrite(
        unmask(bold_noT1gs * Gstd[mask][:, np.newaxis], mask),
        aff,
        "hik_ts_OC_T1c.nii",
        header=head,
    )

    """
    Make medn version of T1 corrected time series
    """
    niwrite(
        Gmu[:, :, :, np.newaxis]
        + unmask((bold_noT1gs + resid) * Gstd[mask][:, np.newaxis], mask),
        aff,
        "dn_ts_OC_T1c.nii",
        header=head,
    )

    """
    Orthogonalize mixing matrix w.r.t. T1-GS
    """
    mmixnogs = mmix.T - np.dot(np.linalg.lstsq(glsig.T, mmix, rcond=None)[0].T, glsig)
    mmixnogs_mu = mmixnogs.mean(-1)
    mmixnogs_std = mmixnogs.std(-1)
    mmixnogs_norm = (mmixnogs - mmixnogs_mu[:, np.newaxis]) / mmixnogs_std[
        :, np.newaxis
    ]
    mmixnogs_norm = np.vstack(
        [np.atleast_2d(np.ones(max(glsig.shape))), glsig, mmixnogs_norm]
    )

    """
    Write T1-GS corrected components and mixing matrix
    """
    sol = np.linalg.lstsq(mmixnogs_norm.T, dat.T, rcond=None)
    niwrite(unmask(sol[0].T[:, 2:], mask), aff, "betas_hik_OC_T1c.nii", header=head)
    np.savetxt("meica_mix_T1c.1D", mmixnogs)

    return glsig
