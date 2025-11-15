# artifact.py
# Utilities for component selection for lib_tedana and related files
# (c) 2025 Prantik Kundu, PhD

from enum import Enum

import numpy as np
import scipy.ndimage as ndimage

from .filter import gradmask, gradmask_sep, localfwhm_afni, spatclust
from .volume import place3d

ARTIFACT_TYPES = Enum(
    "ARTIFACT_TYPES",
    [
        "VESSEL",
        "FOURIER",
    ],
)


def score_vessel_artifact(ted, head, aff):
    # Get data and mask
    psc = ted[..., 0]
    # fs0 = ted[..., 2]
    # zmap = ted[..., 3]
    maskvol = psc != 0
    fwhm = localfwhm_afni(
        psc,
        maskvol,
        header=head,
        affine=aff,
        domask=False,
        rmoutfile=True,
        # prefix="fwhm",
    )
    src_anchor = np.array(maskvol.shape) / 2

    # Cluster parameters err fwhm
    fwhm_min = fwhm.min(-1) == -1
    nbrvox = maskvol.sum()
    clfac = 0.001
    clust_size = int(clfac * nbrvox)
    fwhm_min_clust = spatclust(fwhm_min, thr=0, csize=clust_size) != 0

    fwhm_displ = (fwhm_min_clust).astype(float).copy()  # type: ignore
    for _stride in range(-4, 5):
        fwhm_displ = place3d(
            fwhm_min_clust.astype(int),  # type: ignore
            src_anchor,
            fwhm_displ,
            src_anchor + (0, 0, _stride),
            add=True,
        )
        fwhm_displ = place3d(
            fwhm_min_clust.astype(int),  # type: ignore
            src_anchor,
            fwhm_displ,
            src_anchor + (0, _stride, 0),
            add=True,
        )
        fwhm_displ = place3d(
            fwhm_min_clust.astype(int),  # type: ignore
            src_anchor,
            fwhm_displ,
            src_anchor + (_stride, 0, 0),
            add=True,
        )
    closing_struct = ndimage.generate_binary_structure(3, 2)
    fwhm_displ = ndimage.binary_closing(fwhm_displ, closing_struct).astype(int)
    fwhm_displ = fwhm_displ - fwhm_displ * fwhm_min_clust

    # Separate grads
    gradpsc = gradmask_sep(psc, closingstep=1)
    gradpsc_mask = (gradpsc.max(-1) != 0).astype(int)
    olap = fwhm_displ + 2 * gradpsc_mask
    err_frac = (olap == 3).astype(int).sum() / gradpsc_mask.sum()

    return err_frac


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


def artifact_score_penalty(Kappas, Rhos, ci):
    K = Kappas[ci]  # ctab[ci, 1]
    R = Rhos[ci]  # ctab[ci, 2]
    return R / K**1.25
