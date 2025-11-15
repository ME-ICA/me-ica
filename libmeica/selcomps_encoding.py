# selcompsspectral.py
# Spatial-spectral component selection of in-plane acceleration artifacts and
#   physiological artifacts.
# (c) 2025 Prantik Kundu, PhD

import hashlib
from functools import cached_property
from pathlib import Path
from time import time

import numpy as np
from scipy.stats import stats

from libmeica.utils.selection import (
    andb,
    dice,
    getelbow,
    getelbow2,
    getelbow3,
    getfbounds,
)

from .selcomps_base import SelcompsBase, _u
from .utils.artifact import score_fourier_artifact_count

ACORR_NEG_COUNT_FILE = "acorr_neg_count_file.1D"
ENCODING_FSPACE_FILE = "encoding_fspace.1D"
ENCODING_FSPACE_GT_FILE = "encoding_fspace.1D.gt"

_PRESEL = None


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
    if sample is None and _PRESEL is not None:
        sample = a[_PRESEL]
    elif sample is None:
        sample = a
    _mad = stats.median_abs_deviation(sample)  # type: ignore
    _med = np.median(sample)
    return (a - _med) / _mad


def z1(a, *, sample=None):
    _za = z_(a, sample=sample)
    _za[_za < 1] = 1
    return _za


def zu(a, floor=None, *, sample=None):
    _zu = z_(a, sample=sample) + 1
    _zu[_zu >= 0] = _zu[_zu >= 0] + 1
    _zu[_zu < 0] = np.exp(_zu)[_zu < 0]
    if floor is not None:  # TODO: This is redundant with z1
        _zu[_zu < floor] = floor
    return _zu


def mad_sel(a, sigma=5, sample=None):
    if sample is None and _PRESEL is not None:
        sample = a[_PRESEL]
    elif sample is None:
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
        self.report_elbows()

    @cached_property
    def fourier_artifacts_score(self, force=False):
        # if not force and Path(ACORR_NEG_COUNT_FILE).exists():
        #     err_count_rows = np.loadtxt(ACORR_NEG_COUNT_FILE)
        # else:
        err_count_rows = np.zeros([self.nc, 3])
        err_mag_rows = np.zeros([self.nc, 1])
        for _c in self.nc_:
            err_count_rows[_c], err_mag_rows[_c] = score_fourier_artifact_count(
                self.as_ted(_c), self.header, self.affine
            )
            # np.savetxt(ACORR_NEG_COUNT_FILE, err_count_rows)
        return err_count_rows, err_mag_rows

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
        fmts[14] = "%2.2f"

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
            "Zu:K-R",
            "Zu:r2-s0",
            "Zu:ex",
            "Z1:er",
            "Zu(0.5):R",
            "K-R",
            "r2-s0",
            "ex",
            "er",
            "R",
            "K",
            "Vx",
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

        rej = self.rej

        _n = self.nc_

        K = self.Kappas
        R = self.Rhos
        sr2 = self.countsigFR2
        ss0 = self.countsigFS0

        # import pudb
        #
        # pudb.set_trace()

        _ex, _em = self.fourier_artifacts_score
        ex = _ex.sum(1)
        em = _em.flatten()
        er = z_(_ex.max(1)) / z1(_ex.min(1))

        fmin, fmid, _ = getfbounds(self.Ne)  # type: ignore
        _K_R = K - R
        _sr2_ss0 = sr2 - ss0

        zu_KR = zu(_K_R)
        zu_sr2_ss0 = zu(sr2 - ss0)
        zu_ex = zu(ex)
        zu_er = zu(er)
        zu05_R = zu(R, 0.5)

        denom_max = np.array(
            [
                zu_ex,
                zu_er,
                zu05_R,
            ]
        ).max(0)

        score = Mg(zu_KR, zu_sr2_ss0) / Mg(
            zu_ex,
            zu_er,
            zu05_R,
        )

        knob = 50

        min_K = np.min([getelbow2(K, True), np.average(K, weights=1.0 / self.varex)])
        max_zuex = np.average(zu_ex, weights=R)

        acc_cand = rs[
            andb([score[rs] > 2, K[rs] > min_K, denom_max[rs] < max_zuex]) == 3
        ]

        keep = np.union1d(
            acc_cand[score[acc_cand] > zu_KR[acc_cand]],
            acc_cand[score[acc_cand] >= getelbow(score, True)],
        )

        tent = np.setdiff1d(acc_cand, keep)
        R_thr = np.mean([np.percentile(R[keep], 100 - knob), fmin])
        ex_thr = np.percentile(ex[keep], knob)
        score_thr = np.max(
            [
                np.min(score[keep]),
                np.mean([getelbow(score, True), getelbow2(score, True)]),
            ]
        )

        eject = tent[
            andb([R[tent] > R_thr, ex[tent] > ex_thr, score[tent] < score_thr]) == 3
        ]

        acc = np.setdiff1d(acc_cand, eject)

        rs = np.setdiff1d(rs, acc)
        midk = rs[zu(ex[rs]) > 5]

        ign = np.setdiff1d(rs, midk)

        # Make take of decision process
        markers = np.ones(self.nc) * 9
        markers[acc] = 1
        markers[rej] = 0
        markers[midk] = -1
        markers[ign] = 3
        markers[eject] = 2

        _gt = self.ground_truth()

        self.Fspace = [
            _n,
            markers,
            _gt,
            score,
            zu(_K_R),
            zu(_sr2_ss0),
            zu_ex,
            zu_er,
            zu05_R,
            _K_R,
            _sr2_ss0,
            ex,
            er,
            R,
            K,
            self.varex,
        ]

        self._write_table()

        return acc, rej, midk, ign
