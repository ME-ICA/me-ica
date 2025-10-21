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
        ex = self.fourier_artifacts_score.sum(1)
        er = z_(self.fourier_artifacts_score.max(1)) / z1(
            self.fourier_artifacts_score.min(1)
        )

        fmin, fmid, _ = getfbounds(self.Ne)  # type: ignore
        _K_R = K - R
        _sr2_ss0 = sr2 - ss0

        zu_KR = zu(_K_R)
        zu_sr2_ss0 = zu(sr2 - ss0)
        zu_ex = zu(ex)
        zu_er = zu(er)
        zu05_R = zu(R, 0.5)

        zu_KR_v_ex = zu_KR / zu_ex

        score = Mg(zu_KR, zu_sr2_ss0) / Mg(
            zu_ex,
            zu_er,
            zu05_R,
        )

        # import pudb; pudb.set_trace()

        score_thr = np.mean([5, np.median(zu_KR[K > getelbow2(K, True)])])

        # recognizing that gradient artifact may be kappa-coupled
        #   and Rho may be Kappa coupled, we need to pick a representative
        #   guess for K that discounts off the coupling factors
        #   while still holding the assumption that good high K components
        #   will have basically low zu_ex, even if not as low as the empty
        #   components or pure Rho components
        zu_ex_thr = np.average(zu_ex, weights=zu_KR_v_ex)

        # In high dimensional datasets, 1./self.varex pushes the R threshold
        # tigher but in noisy datastes R might still have high kappa coupling
        # so we pick the more inclusive of the two conditions, the average R
        # consistent with high K values, or the R corresponding to the floor
        R_thresh = np.mean([np.average(R, weights=zu_KR_v_ex), fmin, getelbow(R, True)])

        acc_high = _n[andb([score > score_thr, zu_ex_thr > zu_ex, R < R_thresh]) == 3]

        # Make an ignore set that will be whittled down
        ign = np.setdiff1d(rs, acc_high)

        # Get rid of high-variance clear artifacts
        if len(acc_high) >= 4:
            _sample = acc_high
        else:
            _sample = rs
        zu_ex_art_thr = np.mean([np.average(zu_ex, weights=R), 5])
        midk_cand_1 = ign[
            mad_sel(self.varex[ign], sigma=4, sample=self.varex[_sample])
        ]  # High variance
        midk_cand_2 = _n[zu_ex >= zu_ex_art_thr]  # Clear artifact
        midk = np.intersect1d(midk_cand_1, midk_cand_2)

        # Whittle
        rs = np.setdiff1d(ign, midk)

        # Keep subtle TE-dependent components
        if len(acc_high) >= 4:
            _sample = acc_high
        else:
            _sample = rs
        zu_ex_thr_2 = np.mean([np.average(zu_ex, weights=K), np.max(zu_ex[_sample])])
        # zu_er_thr_2 = np.mean([np.average(zu_er, weights=K), np.max(zu_er[acc_high])])
        R_thresh_2 = np.mean([fmid, np.max(R[_sample])])

        # The _2 tresholds are relaxed compared to the acc_high pass
        print(f"""
        >score_thr={score_thr}
        <zu_ex_thr={zu_ex_thr}
        <R_thresh={R_thresh}
        <<zu_ex_art_thr={zu_ex_art_thr}
        <zU_ex_thr_2={zu_ex_thr_2}
        <R_thresh_2={R_thresh_2}
        """)

        # Basic inclusion
        keep_1 = rs[
            andb(
                [
                    # Guessing location of 'ideal' elbow
                    K[rs] > getelbow2(K, True),
                    score[rs] > 3,
                    zu_er[rs] < 1,
                    zu_ex[rs] <= zu_ex_thr_2,
                    R[rs] <= R_thresh_2,
                ]
            )
            == 5
        ]
        # Safety condition to high-res / 7T where er drives a bit too much exlusion
        keep_2 = rs[
            andb(
                [
                    K[rs] > getelbow2(K, True),
                    score[rs] > 4,
                    zu_ex[rs] < zu_ex_thr_2,
                    R[rs] < R_thresh_2,
                ]
            )
            == 4
        ]
        keep = np.union1d(keep_1, keep_2)

        # Final
        acc = np.union1d(acc_high, keep)  # to go into Hi-K
        # Extra harmless variance
        ign = np.setdiff1d(ign, np.union1d(np.union1d(keep, midk), acc_high))

        # Make take of decision process
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
