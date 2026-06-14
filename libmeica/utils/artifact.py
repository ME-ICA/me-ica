# artifact.py
# Utilities for component selection for lib_tedana and related files
# (c) 2025 Prantik Kundu, PhD

import numpy as np

from .filter import gradmask, localfwhm_afni


def score_fourier_artifact_count(ted, head, aff, prefix=""):
    psc = ted[..., 0]
    maskvol = psc != 0
    fwhm = localfwhm_afni(
        psc,
        maskvol,
        header=head,
        affine=aff,
        domask=False,
        prefix=prefix,
        # rmoutfile=True,
    )
    fwhm_sel = fwhm == -1
    fwhm_count = fwhm_sel.sum(0).sum(0).sum(0)

    fwhm_max = fwhm.max(-1)
    fwhm_min = fwhm.min(-1)
    psc[psc < 0] = 0

    ign = fwhm_min == -1
    fwhm_psc = np.average(fwhm_max[~ign], weights=psc[~ign] ** 2.0)
    fwhm_fr2 = np.average(fwhm_max[~ign], weights=ted[..., 1][~ign] ** 2.0)
    fwhm_fs0 = np.average(fwhm_max[~ign], weights=ted[..., 2][~ign] ** 2.0)

    fwhm_sig = np.array([fwhm_psc, fwhm_fr2, fwhm_fs0])
    return fwhm_count, fwhm_sig


def score_fourier_artifact(ted, head, aff):
    psc = ted[..., 0]
    maskvol = psc != 0
    fwhm = localfwhm_afni(
        psc,
        maskvol,
        header=head,
        affine=aff,
        domask=False,
        rmoutfile=True,
    )
    gradpsc = gradmask(psc, norm=False)
    fwhm_sel = fwhm == -1
    err_frac = (fwhm_sel * np.abs(gradpsc)).sum()
    return err_frac
