# selcompsspectral.py
# Spatial-spectral component selection of in-plane acceleration artifacts and
#   physiological artifacts.
# (c) 2025 Prantik Kundu, PhD

from pathlib import Path
from typing import List, Tuple

import numpy as np

from .utils.filter import niwrite
from .utils.selection import andb, getelbow, getelbow2, getfbounds
from .utils.volume import unmask

_u = np.union1d


class SelcompsBase:
    """
    This class takes as input a seldict from fitmodels and its
    component table and performs spectral analysis of spatial components
    in concert with Kappa and Rho values in order to segregate in-plane
    acceleration artifacts and to some extent physiological artifact
    from clearly neural components.
    Neural components here are defined as:
        - TE-dependent (High-Kappa)
        - Cauchy or 1/f -like,
            - i.e. sharply higher amplitude at lower frequency
            - and spectrally sparse otherwise
        - Spectrally isotropic, e.g. no preferred spatial frequency or angle.
    See provisional patent filing for more discussion.
    """

    def __init__(
        self, *, seldict, Ne, mask, t2s, header, affine, midk_fac=1, curtail=True
    ):
        self.seldict = seldict
        # self.ctab = ctab
        self.midk_fac = midk_fac
        self.Ne = Ne
        self.header = header
        self.affine = affine
        # mask and t2sG will be pulled from a T2s object to be implemented
        self.mask = mask
        self.t2s = t2s  # make sure t2s is smaller than mask, which is same size as t2sG
        self._rej = None
        self.curtail = False
        if getelbow(self.Kappas) * 3 < self.nc and curtail:
            self.curtail = True
        self.mmix_id = seldict["mmix_id"]

        # Pull in maps
        fmin, _, _ = getfbounds(self.Ne)
        countFS0 = self.seldict["F_S0_maps"] > fmin
        countFR2 = self.seldict["F_R2_maps"] > fmin
        countZ = self.seldict["Z_maps"] > 1.95
        self.countFS0 = countFS0.sum(0)
        self.countFR2 = countFR2.sum(0)
        self.countZ = countZ.sum(0)
        self.FS0lclsz = self.seldict["F_S0_lclsz"]
        self.FR2lclsz = self.seldict["F_R2_lclsz"]
        self.Zlclsz = self.seldict["Z_lclsz"]
        self.countsigFS0 = self.seldict["F_S0_clmaps"].sum(0)
        self.countsigFR2 = self.seldict["F_R2_clmaps"].sum(0)
        self.countsigZ = self.seldict["Z_clmaps"].sum(0)
        self.corrT1mi = seldict["corrT1mi"]
        self.corrT2s = seldict["corrT2s"]

    # Here are some utility functions
    @property
    def Kappas(self):
        return self.seldict["Kappas"]
        # return self.ctab[:, 1]

    @property
    def Rhos(self):
        return self.seldict["Rhos"]
        # return self.ctab[:, 2]

    @property
    def varex(self):
        return self.seldict["varex"]

    @property
    def mmix(self):
        return self.seldict["mmix_new"]

    @property
    def rej(self):
        """
        Calculate the rejected set, cache, and return the array
        """
        if self._rej is None:
            self._rej = self.nc_[
                andb([self.Rhos > self.Kappas, self.countsigFS0 > self.countsigFR2]) > 0
            ]
        return self._rej

    @property
    def nc(self):
        """
        Calculate the total dimensionality of the data
        """
        return self.seldict["Kappas"].shape[0]

    @property
    def nc_(self):
        """
        Return the list of all indicies
        """
        return np.arange(self.nc)

    @property
    def ncl(self):
        """
        Calculate and return the indices of components that are not rejected
        or at the ignored tail of the Kappa spectrum
        """
        if self.curtail:
            ign = self.tail_to_ign(self.nc_, rho=False)
            return np.setdiff1d(self.nc_, _u(self.rej, ign))
        else:
            return np.setdiff1d(self.nc_, self.rej)

    @property
    def n4c(self):
        """Return the number of components to classify post-rejection
        and tail ignore"""
        return len(self.ncl)

    def as_ted(self, i):
        """
        Return a 4D array with standard briks as outputted in cc***.nii.gz
        files
        """
        return np.array(
            [
                unmask(self.seldict["PSC"][:, i], self.mask),
                unmask(self.seldict["F_R2_maps"][:, i], self.t2s != 0),
                unmask(self.seldict["F_S0_maps"][:, i], self.t2s != 0),
                unmask(self.seldict["Z_maps"][:, i], self.mask),
            ]
        ).transpose(1, 2, 3, 0)

    def write_ted(self, i, fn=None):
        if fn is None:
            outpath = Path("/tmp") / ("ted_%03i.nii" % i)
        else:
            outpath = Path(fn)
        niwrite(self.as_ted(i), self.affine, outpath, self.header)
        return outpath

    def midk_detected(self, two_ks_mat) -> bool:
        """
        Detect if there is a consensus decision function for all comparisons
        of each component against all the others. If there's a small number of
        cut-points that all the components show in comparison with all others,
        then there's a midk group.
        """

        raise NotImplementedError

        # nc = self.n4c
        # fac = self.midk_fac
        # print(f"Midk fac is {nc/fac}")
        # if len(np.unique([getelbow3(two_ks_mat[i]) for i in range(nc)])) < nc / fac:
        #     # e.g. Do we have a small number of detected thresholds?
        #     return True
        # else:
        #     return False

    def sel_midk(self) -> Tuple[List, List, List]:
        """
        This function takes all components not captured as rej and groups
        them into accepted, midk, and an ignore group defined on not fitting
        either of those. The returned group components account for all
        components.
        """
        raise NotImplementedError

        # return acc, midk, ign_marg

    def tail_to_ign(self, ingrp, rho=False):
        """
        Computes Kappa and optionally Rho thresholds conservatively and
        reports back the indices within ingrp that are below the tresholds
        """
        # Pick the Rho and Kappa thr representing most conservative estimates
        Rhos_thr = np.mean(
            [getelbow(self.Rhos, True), getelbow2(self.Rhos, True)]
            + list(getfbounds(4)[:2])
        )
        Kappas_thr = np.mean([getelbow(self.Kappas, True)] + list(getfbounds(4)[:2]))
        if rho:
            ign_rho = ingrp[self.Rhos[ingrp] > Rhos_thr]
        else:
            ign_rho = []
        to_clf = np.setdiff1d(ingrp, ign_rho)
        ign_kappa = to_clf[self.Kappas[to_clf] < Kappas_thr]
        return _u(ign_kappa, ign_rho)

    def good_guess(self):
        good_guess_ted = self.nc_[
            andb(
                [
                    self.Kappas > getelbow(self.Kappas, True),
                    self.Rhos < getelbow(self.Rhos, True),
                ]
            )
            == 2
        ]
        return good_guess_ted

    def report_elbows(self):
        print(
            getelbow(self.Kappas, True),
            getelbow2(self.Kappas, True),
            getelbow(self.Rhos, True),
            getelbow2(self.Rhos, True),
            getfbounds(self.Ne),
        )

    def fit(self):
        """
        Computes the spectral Kolomogorov-Smirnov distances between
        component spatial spectra to distinguish spectrally sparse, isotropic,
        low-frequency components in order to group neural components separately
        from physiological components (which have preferred angles), in-plane
        acceleration artifacts (which have high-frequencies), and
        statistical noise (which are not sparse)
        """
        midk = self.sel_midk()

        acc = np.setdiff1d(self.ncl, _u(midk, self.rej))

        ign = np.setdiff1d(self.nc_, _u(_u(acc, midk), self.rej))

        return acc, self.rej, midk, ign
