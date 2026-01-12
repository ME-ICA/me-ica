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

from numpy import ptp, min, percentile

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

pct = lambda p,v: percentile(p,v,method='linear')

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


def zu(a, floor=None, *, bias=1, sample=None):
    # bias creates a key dynamic range behavior
    # it keeps data < 1 MAD from median in equal consieration
    # which is fair in noisy data
    # and allows zu(x) and zu(-x) to create clamps/grips
    _zu = z_(a, sample=sample) + bias
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
        err_mag_rows = np.zeros([self.nc, 3])
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

    def _write_table(self, header=""):
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

        full_header = header + '\n' + "".join(label_string)

        np.savetxt(
            ENCODING_FSPACE_FILE,
            np.array(self.Fspace).T,
            fmt=fmts,
            header=full_header,
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

        # Simple optimization to discover ex_thr
        def r2_opt(_sr2, _ex, *, cost=1, slack=0.8):
            if cost == 1:
                r2_cost = zu(_sr2) - z_(_ex)
            else:
                r2_cost = zu(_sr2) - zu(_ex)
            ex_costs = []
            for _exv in sorted(_ex):
                ex_costs.append(r2_cost[_ex <= _exv].sum())
            ex_thr = np.sort(_ex)[np.arange(len(_ex))[ex_costs / np.max(ex_costs) >= slack].max()]
            ex_thr_ns = np.sort(_ex)[np.arange(len(_ex))[ex_costs / np.max(ex_costs) >= 1].max()]
            return ex_thr, ex_thr_ns, ex_costs 

        # Set K thresholds
        K_thrs = np.array(
            [
                getelbow(K, True),
                getelbow2(K, True),
                np.average(K, weights=1.0 / self.varex),
                fmin,
                fmid,
            ]
        )
        K_thr_safe, K_thr_agg  = np.percentile(K_thrs, [25,50], method='linear')
        
        # Determine dataset conditions
        high_grappa = False
        high_rho = False
        if (K.max() < fmax and len(rs)/len(self.nc_) > 0.75):
            high_grappa = True
            # This if-condition shows the dataset to have low max K, yet K>R almost everywhere
            #  suggesting a dataset affectd by heavy GRAPPA/T2* weighted artifact
            #  so we raise the Kappa thresholds
            K_thr_safe, K_thr_agg  = np.percentile(K_thrs, [50,75], method='linear')

        # High Rho dataset
        if getelbow2(R, True) > getelbow2(K, True):  
            high_rho = True
        R_thr = np.min([np.mean([getelbow(R, True), getelbow2(R,True)]), 1.5*fmin])

        # Set optimization set
        if high_grappa or high_rho:
            # Determine ex thr from enriched set
            rs_opt = self.nc_[score>getelbow2(score,True)]
            rs_opt = np.intersect1d(rs_opt, rs)
        else:
            rs_opt = rs

        # Define aggressive and conservative ex thresholds
        ex_thr_1, ex_thr_1ns, _ex_thr_costs_1 = r2_opt(sr2[rs_opt], ex[rs_opt], cost=1)
        ex_thr_2, ex_thr_2ns, _ex_thr_costs_2 = r2_opt(sr2[rs_opt], ex[rs_opt], cost=2)
        ex_thr_agg = np.min([ex_thr_1, ex_thr_2])
        if high_grappa or high_rho:
            ex_thr_agg = np.min([ex_thr_1ns, ex_thr_2ns])

        # Select calibration set
        acc_sel_a = andb(
                [
                    K[rs] >= K_thr_agg,
                    ex[rs] < ex_thr_agg,
                    R[rs] <= R_thr,
                ]
            )
        acc_a = rs[acc_sel_a==3]
        eject_a_ss0 = acc_a[rankvec(ss0[acc_a]) - rankvec(K[acc_a]) >= np.max([6,len(acc_a)/2])]
        eject_a_R = acc_a[rankvec(R[acc_a]) - rankvec(K[acc_a]) >= np.max([6,len(acc_a)/2])]
        eject_a_er = acc_a[rankvec(er[acc_a]) - rankvec(K[acc_a]) >= np.max([6,len(acc_a)/2])]
        eject_a = np.union1d(eject_a_R, eject_a_ss0)
        eject_a = np.union1d(eject_a, eject_a_er)
        acc_a = np.setdiff1d(acc_a, eject_a)

        # Nonparametric thresholds
        sr2_ex = sr2/ex
        sr2_ss0 = sr2-ss0
        sr2_ex_thr = sr2_ex[acc_a].min() * 0.2
        sr2_ss0_thr = sr2_ss0[acc_a].min() * 0.2
        sr2_thr = np.percentile(sr2[acc_a], 20) * 0.2

        acc_sel_b = andb(
                [
                sr2[rs] > sr2_thr,
                sr2_ex[rs] > sr2_ex_thr,
                sr2_ss0[rs] > sr2_ss0_thr,
                K[rs] >= K_thr_safe,
            ]
        )
        acc_b = rs[acc_sel_b == 4]

        # Clean up based on R scores in acc_b while protecting robust R2* task effects
        z_sr2_ss0_a = ((sr2-ss0)[_n]-(sr2-ss0)[acc_a].mean())/(sr2-ss0)[acc_a].std()
        z_R_a = (R[_n]-R[acc_a].mean())/R[acc_a].std()
        eject_R = acc_b[andb([(z_R_a - z_sr2_ss0_a)[acc_b] > 0, z_R_a[acc_b] > 0])==2]

        # Remove based on high relative ex compared to acc_a
        z_ex_a = (ex[_n]-ex[acc_a].mean())/ex[acc_a].std()
        z_sr2_ex_a = (sr2_ex[_n]-sr2_ex[acc_a].mean())/sr2_ex[acc_a].std()
        eject_ex = acc_b[andb([(z_ex_a - z_sr2_ss0_a)[acc_b] > 6.5, z_ex_a[acc_b] > 4])==2]

        # High anisotropy and large magnitude of ACF indicates a low-freq
        #  gradient field artifact. Extreme anisotropy can happen at middle of 
        #  K-spectrum, so no other conditions. Moderate happens at tail so can
        #  control using K_thr
        eject_field = np.intersect1d(
            rs[z_(np.max(_em, axis=1)[rs]) > 9],
            rs[z_((np.max(_em, axis=1) / np.min(_em, axis=1))[rs]) > 9],
        )
        eject_field_rank = acc_b[rankvec(_em.max(axis=1)[acc_b]) - rankvec(K[acc_b]) > len(acc_a)]
        eject_field = np.union1d(eject_field, eject_field_rank)
        
        # Remove large included artifacts after tail inclusion
        eject_score_thr = np.mean([getelbow(score,True), np.percentile(score[acc_a], 33)])
        high_varex_thr = np.mean(self.varex[acc_a])
        eject_quality = acc_b[andb([score[acc_b]<eject_score_thr, R[acc_b]>R_thr, self.varex[acc_b]>high_varex_thr])>=2]
        
        # Combine ejects
        eject_R = np.setdiff1d(eject_R, acc_a)  # R not allowed to go after acc_a
        eject = np.union1d(eject_field, eject_R)  # field will be protected by quality requirements
        eject = np.intersect1d(eject, eject_quality)
        eject_ex = np.setdiff1d(eject_ex, acc_a)
        eject = np.union1d(eject, eject_ex)

        # Final acceptance
        acc = np.setdiff1d(acc_b, eject) 

        min_score = score[acc].min()
        print("Minimum score is", min_score)

        rs = np.setdiff1d(rs, acc)
        midk = rs[zu_ex[rs] > 5]

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

        threshold_string = f"K_thrs: {K_thr_safe} {K_thr_agg}, score thrs: {eject_score_thr}, ex thrs: {ex_thr_agg}, n_eject:{len(eject)}"
        self._write_table(header=threshold_string)

        return acc, rej, midk, ign
