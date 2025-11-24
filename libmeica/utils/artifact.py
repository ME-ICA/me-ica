# artifact.py
# Utilities for component selection for lib_tedana and related files
# (c) 2025 Prantik Kundu, PhD

import numpy as np

from .filter import gradmask, localfwhm_afni


def score_fourier_artifact_count(ted, head, aff):
    psc = ted[..., 0]
    maskvol = psc != 0
    fwhm = localfwhm_afni(
        psc,
        maskvol,
        header=head,
        affine=aff,
        domask=False,
        rmoutfile=False,
    )
    fwhm_sel = fwhm == -1
    fwhm_count = fwhm_sel.sum(0).sum(0).sum(0)

    sel_any = fwhm_sel.sum(3)

    mag = np.abs(psc)
    sig = np.abs(ted[..., 3])
    sig[sig < 1] = 1

    frac_power = (sel_any * mag * sig).sum() / mag.sum()

    # fwhm_mean = np.average(np.abs(psc) ** 2, weights=sel_any) ** 1.0 / 2

    return fwhm_count, frac_power


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
