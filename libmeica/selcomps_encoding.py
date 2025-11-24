# selcompsspectral.py
# Spatial-spectral component selection of in-plane acceleration artifacts and
#   physiological artifacts.
# (c) 2025 Prantik Kundu, PhD

import hashlib
from functools import cached_property
from pathlib import Path
from sys import argv
from time import time


import numpy as np
from scipy.stats import stats

from libmeica.utils.selection import (
    andb,
    dice,
    getelbow,
    getelbow2,
    getfbounds,
)

from .selcomps_base import SelcompsBase
from .utils.artifact import score_fourier_artifact_count
from .utils.filter import rankvec

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

        _ex, _em = self.fourier_artifacts_score
        ex = _ex.sum(1)
        er = z_(_ex.max(1)) / z1(_ex.min(1))

        if "DEBUG" in argv:
            import pudb

            pudb.set_trace()

        fmin, fmid, fmax = getfbounds(self.Ne)  # type: ignore
        _K_R = K - R
        # _R_K = R - K
        _sr2_ss0 = sr2 - ss0

        zu_KR = zu(_K_R)

        zu_sr2_ss0 = zu(sr2 - ss0)
        zu_ex = zu(ex)
        zu_er = zu(er)
        zu05_R = zu(R, 0.5)

        score = Mg(zu_KR, zu_sr2_ss0) / Mg(
            zu_ex,
            zu_er,
            zu05_R,
        )

        K_thr = np.mean(
            [
                getelbow(K, True),
                getelbow2(K, True),
                np.average(K, weights=1.0 / self.varex),
            ]
        )

        R_thrs = [
            # getelbow(R, True),
            np.average(R, weights=1.0 / self.varex),
            np.average(R, weights=1.0 / K),
            np.average(R[rej], weights=1.0 / R[rej]),
        ]
        R_thr = np.mean(R_thrs)

        # Make a good guess of a set of decent components
        acc_ini = []

        # Guess how many good components there are
        est_n_good = (
            andb([K[rs] > getelbow2(K, True), zu_ex[rs] < 3, R[rs] < R_thr]) == 3
        ).sum()

        # Make an initial guess of a medium score
        score_start = np.max([np.percentile(score, 95) / 10, 7])
        score_ini_thr = score_start
        for _s in range(100):
            score_ini_thr = score_start - _s / 10
            acc_ini = rs[score[rs] > score_ini_thr]
            # acc_ini = acc_ini[R[acc_ini] < R_thr]
            if len(acc_ini) > est_n_good / 2:
                break

        # Check if a initial guess was bad and place stopgap
        if score_ini_thr < 1 or len(acc_ini) <= 4:
            score_ini_thr = 3  # getelbow(score, True)
            acc_ini = score > score_ini_thr
            print(
                "WARNING Couldn't find a stable initial guess for component selection, results likely wrong. "
            )

        # Estimate the ex thr (GRAPPA artifact) from the iniital guess
        zu_nex = zu(-ex)
        zu_ner = zu(-er)
        zu_nex_thr = np.min(
            [2, zu_nex[acc_ini].min()]
        )  # Note this does not scale like zu_ex, bigger is better!
        zu_ner_thr = np.min([zu_nex_thr / 2, zu_ner[acc_ini].min() / 2])

        # Final selection from a group of rational votes
        sel_sum = andb(
            [
                zu_KR > 2,
                zu(K) > 2,
                zu_sr2_ss0 > 2,
                K >= K_thr,
                # Dynamic bounds
                score > 2,
                zu_nex >= zu_nex_thr,
                zu_ner > zu_ner_thr,
                # Hard lower bounds
                zu_nex >= np.min([zu_nex_thr, 1]),
                zu_ner > 0,
            ]
        )

        if score_ini_thr < 5 and est_n_good > 20:
            acc_wR = np.intersect1d(_n[sel_sum >= 7], rs)
            print(
                "Warning! Found high mixed artifact effect or nigh neural load and low artifact load, falling back to conservative component selection."
            )
        else:
            acc_wR = np.intersect1d(_n[sel_sum == 9], rs)

        # Heuristic exclusion of "high" R comps to ignore values on their MAD
        eject = acc_wR[
            andb(
                [
                    rankvec(R[acc_wR]) - rankvec(zu_KR[acc_wR]) > est_n_good / 2,
                    R[acc_wR] > R_thr,
                    zu(R[acc_wR]) > 6,
                ]
            )
            == 3
        ]

        acc = np.setdiff1d(acc_wR, eject)

        min_score = score[acc].min()
        print("Minimum score is", min_score)

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
