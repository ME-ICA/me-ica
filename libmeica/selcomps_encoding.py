# selcompsspectral.py
# Spatial-spectral component selection of in-plane acceleration artifacts and
#   physiological artifacts.
# (c) 2023 Prantik Kundu, PhD
# Ceretype Neuromedicine Inc.

import hashlib
from functools import cached_property
from pathlib import Path
from time import time

import numpy as np
from scipy.stats import stats

from libmeica.utils.selection import andb, dice, getelbow2

from .selcomps_base import SelcompsBase, _u
from .utils.artifact import score_fourier_artifact_count

ACORR_NEG_COUNT_FILE = "acorr_neg_count_file.1D"
ENCODING_FSPACE_FILE = "encoding_fspace.1D"
ENCODING_FSPACE_GT_FILE = "encoding_fspace.1D.gt"


def zp(a):
    if np.median(a) == 0:
        # breakpoint()
        a = a.copy()
        _end = np.min(a[a > 0])
        _start = np.max(a[a <= 0])
        _cnt = int(np.sum(a == 0))
        a[a == 0] = np.linspace(_start, _end, _cnt)
    return a


def z_(a, *, sample=None):
    a = zp(a)
    if sample is None:
        sample = a
    _mad = stats.median_abs_deviation(sample)  # type: ignore
    _med = np.median(sample)
    return (a - _med) / _mad


def z1(a, *, sample=None):
    _za = z_(a, sample=sample)
    _za[_za < 1] = 1
    return _za


def zu(a, floor=0.01, *, sample=None):
    _zu = z_(a, sample=sample)
    _zu[_zu < -1] = np.exp(_zu)[_zu < -1]
    _zu[_zu < floor] = floor
    return _zu


def mad_sel(a, sigma=5, sample=None):
    if sample is None:
        sample = a
    _mad = stats.median_abs_deviation(sample)  # type: ignore
    _med = np.median(sample)
    a_sel = a > (_med + sigma * _mad)
    return a_sel


def generate_md5_from_text(file_path):
    hash_md5 = hashlib.md5()
    with open(file_path, "r") as f:
        # Read the entire file content as a single string
        file_content = f.read()
        # Encode the string to bytes and update the hash object
        hash_md5.update(file_content.encode("utf-8"))
    return hash_md5.hexdigest()


def Mg(*a):
    # breakpoint()
    _Mg = np.power(np.prod(a, axis=0), 1.0 / len(a))
    _Mg[np.isnan(_Mg)] = 0
    return _Mg


class SelcompsEncoding(SelcompsBase):
    def __init__(self, **kwargs):
        self.ts = time()
        super().__init__(**kwargs)

    @cached_property
    def fourier_artifacts_score(self, force=False):
        if not force and Path(ACORR_NEG_COUNT_FILE).exists():
            err_count_rows = np.loadtxt(ACORR_NEG_COUNT_FILE)
        else:
            err_count_rows = np.zeros([self.nc, 3])
            for _c in self.nc_:
                err_count_rows[_c] = score_fourier_artifact_count(
                    self.as_ted(_c), self.header, self.affine
                )
            np.savetxt(ACORR_NEG_COUNT_FILE, err_count_rows)
        return err_count_rows

    @cached_property
    def mmix_hash(self):
        tmp_mmix = f"_mmix_tmp_{self.ts}.1D"
        np.savetxt(tmp_mmix, self.mmix)
        return generate_md5_from_text(tmp_mmix)

    def _write_table(self):
        fmts = ["%5.1f" for _ in range(len(self.Fspace))]
        fmts[:3] = ["%5i"] * 3
        fmts[3:6] = ["%9.1f"] * 3
        fmts[9] = "%8i"
        fmts[10:12] = ["%7i"] * 2

        spacing = [6] * len(self.Fspace)
        spacing[0] = 5
        spacing[:3] = [5] * 3
        spacing[3:6] = [10] * 3
        spacing[9] = 8
        spacing[10:12] = [7] * 2

        labels = [
            "#",
            "Acc",
            "GT",
            "Scr=",
            "Z:K-R",
            "Z:r2-s0",
            "/Zu:ex",
            "Zu:er",
            "zu:R",
            "K-R",
            "r2-s0",
            "ex",
            "er",
            "R",
        ]
        label_string = [
            ("{:>%s}" % spacing[_i]).format(_l) for _i, _l in enumerate(labels)
        ]

        hash_footer = f"{self.mmix_hash} // {self.ts}"

        np.savetxt(
            ENCODING_FSPACE_FILE,
            np.array(self.Fspace).T,
            fmt=fmts,
            header="".join(label_string),
            footer=hash_footer,
        )

    def ground_truth(self, gt_fspace_fn=ENCODING_FSPACE_GT_FILE, *, clear=False):
        # breakpoint()
        new_gt = True
        if Path(gt_fspace_fn).exists() and not clear:
            try:
                hashfooter = open(gt_fspace_fn).readlines()[-1]
                _hash = hashfooter.split("//")[0].split("#")[1].strip()
                if _hash == self.mmix_hash:
                    new_gt = False
            except Exception:
                pass
        if not new_gt:
            gt_fspace = np.loadtxt(gt_fspace_fn)
            return gt_fspace[:, 2]
        if new_gt:
            return np.ones(self.nc) * 9

    def ground_truth_test(self, gt_fspace_fn=ENCODING_FSPACE_GT_FILE):
        self.fit
        cl = self.Fspace[1]
        gt = self.ground_truth(gt_fspace_fn)
        _res = {}
        valid_labels = np.setdiff1d(np.unique(gt), [9])  # type: ignore
        for _lab in valid_labels:
            # breakpoint()
            lab_dice = dice(gt == _lab, cl == _lab)
            _res[_lab] = lab_dice
        if _res != {}:
            print("Accuracy test results:", _res)
        return _res

    @cached_property
    def fit(self):
        # Initial removal of rej
        rs = self.nc_
        rs = np.setdiff1d(rs, self.rej)

        # After getting rid of mids, git rid of extreme rhos
        rs_log_rho = np.log(self.Rhos[rs])
        rs_sel = mad_sel(rs_log_rho)
        rej = _u(self.rej, rs[rs_sel])
        rs = np.setdiff1d(rs, rej)

        _n = self.nc_

        K = self.Kappas
        R = self.Rhos
        sr2 = self.countsigFR2
        ss0 = self.countsigFS0
        ex = self.fourier_artifacts_score.sum(1)
        # ezx = self.fourier_artifacts_score_z.max(1)
        er = z_(self.fourier_artifacts_score.max(1)) / z1(
            self.fourier_artifacts_score.min(1)
        )

        score = Mg(zu(K - R), zu(sr2 - ss0)) / Mg(zu(ex, 0.5), z1(er), zu(R, 0.5))

        subelbow = _n[K <= getelbow2(K, True)]
        acc_high = np.setdiff1d(_n[score > 5], np.union1d(rej, subelbow))
        ign = np.setdiff1d(rs, acc_high)

        # Scavenge the ign
        ign = np.setdiff1d(rs, acc_high)
        midk = ign[mad_sel(self.varex[ign])]
        rs = np.setdiff1d(ign, midk)
        keep = rs[
            andb(
                [
                    z1(R[rs]) == 1,
                    z1(ex[rs]) == 1,
                    er[rs] <= 1,
                    K[rs] > getelbow2(K[rs], True),
                ]
            )
            == 4
        ]
        acc = np.union1d(acc_high, keep)
        ign = np.setdiff1d(ign, np.union1d(keep, midk))

        markers = np.ones(self.nc) * 9
        markers[acc] = 1
        markers[rej] = 0
        markers[midk] = -1
        markers[keep] = 2
        markers[ign] = 3

        _gt = self.ground_truth()

        self.Fspace = [
            _n,
            markers,
            _gt,
            score,
            z_(K - R),
            z_(sr2 - ss0),
            zu(ex),
            zu(er),
            zu(R),
            K - R,
            sr2 - ss0,
            ex,
            er,
            R,
        ]

        self._write_table()

        return acc, rej, midk, ign
