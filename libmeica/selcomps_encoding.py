# selcompsspectral.py
# Spatial-spectral component selection of in-plane acceleration artifacts and
#   physiological artifacts.
# (c) 2025 Prantik Kundu, PhD

from datetime import datetime
from functools import cached_property

import numpy as np

from .fit_encoding import fit_encoding
from .selcomps_base import SelcompsBase
from .utils.artifact import score_fourier_artifact_count
from .utils.encoding_table import EncodingTable
from .utils.selection import mmix_hash, z1, z_


class SelcompsEncoding(SelcompsBase):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.report_elbows()

    @cached_property
    def fourier_artifacts_score(self, force=False):
        # if not force and Path(ACORR_NEG_COUNT_FILE).exists():
        #     err_count_rows = np.loadtxt(ACORR_NEG_COUNT_FILE)
        # else:
        err_count_rows = np.zeros([self.nc, 3])
        err_mag_rows = np.zeros([self.nc, 3])
        print("Evaluating GRAPPA Artifact", end="", flush=True)
        for _c in self.nc_:
            err_count_rows[_c], err_mag_rows[_c] = score_fourier_artifact_count(
                self.as_ted(_c), self.header, self.affine
            )
            print(".", end=".", flush=True)
            # np.savetxt(ACORR_NEG_COUNT_FILE, err_count_rows)
        return err_count_rows, err_mag_rows

    @cached_property
    def fit(self):
        # Initial removal of rej
        fs = EncodingTable(Ne=self.Ne, nc=self.nc)
        fs.N = rs = self.nc_
        rs = np.setdiff1d(rs, self.rej)
        rej = self.rej

        # Initialize feature space
        fs.C[rej] = 0

        # Populate with data
        fs.K = self.Kappas
        fs.R = self.Rhos
        fs.sr2 = self.countsigFR2
        fs.ss0 = self.countsigFS0
        fs.sz = self.countsigZ
        fs.tr2 = self.countFR2
        fs.ts0 = self.countFS0
        fs.tz = self.countZ
        fs.lr2 = self.FR2lclsz
        fs.ls0 = self.FS0lclsz
        fs.lz = self.Zlclsz
        fs.V = self.varex
        _ex, _em = self.fourier_artifacts_score
        fs.ex = _ex.sum(1)
        fs.er = z_(_ex.max(1)) / z1(_ex.min(1))
        fs.emx = np.max(_em, axis=1)
        fs.emr = np.max(_em, axis=1) / np.min(_em, axis=1)
        fs.cT1 = self.corrT1mi
        fs.cT2s = self.corrT2s

        fs.save("encoding_fspace.1D")
        etfn = f"encoding_fspace_{self.mmix_id}.1D"
        fs.save(etfn)

        # Get classes
        e = fit_encoding(etfn)

        acc = self.nc_[e.C == 1]
        midk = self.nc_[e.C == -1]
        ign = self.nc_[e.C == 3]
        rej = self.nc_[e.C == 0]

        if np.setdiff1d(np.unique(e.G), [3]).size == 0:
            e.G = e.C

        e.save(etfn)
        e.save("encoding_fspace.1D")

        return acc, rej, midk, ign
