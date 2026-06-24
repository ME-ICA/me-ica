# selcompsspectral.py
# Spatial-spectral component selection of in-plane acceleration artifacts and
#   physiological artifacts.
# (c) 2025 Prantik Kundu, PhD

import os
from sys import argv

import numpy as np

from libmeica.utils.selection import Mg, andb, getelbow, getelbow2, getfbounds, z_, zu

from .utils.encoding_table import EncodingTable
from .utils.filter import rankvec

ENCODING_FSPACE_FILE = "_encoding_fspace.1D"


def fit_encoding(fspace_fn=ENCODING_FSPACE_FILE):
    if "DEBUG" in argv:
        import pudb

        pudb.set_trace()

    e = EncodingTable(Ne=0)
    e.load(fspace_fn)

    N = e.N.astype(int)
    K = e.K
    R = e.R
    C = e.C.astype(int)
    ex = e.ex
    er = e.er
    varex = e.V
    sr2 = e.sr2
    ss0 = e.ss0

    nc = int(e.N.max()) + 1
    rej = N[C == 0]
    rs = N[C != 0].astype(int)
    nc = N.max() + 1

    fmin, fmid, fmax = getfbounds(e.Ne)  # type: ignore
    e.K_R = _K_R = K - R
    e.sr2_ss0 = _sr2_ss0 = sr2 - ss0

    e.zu_KR = zu_KR = zu(_K_R)
    e.zu_sr2_ss0 = zu_sr2_ss0 = zu(_sr2_ss0)
    e.zu_ex = zu_ex = zu(ex)
    e.zu_er = zu_er = zu(er)
    e.zu05_R = zu05_R = zu(R, 0.5)

    e.zu_cT1 = zu(e.cT1)
    e.zu_cT2s = zu(-e.cT2s)

    e.score = score = Mg(zu_KR, zu_sr2_ss0) / Mg(
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
        ex_thr = np.sort(_ex)[
            np.arange(len(_ex))[ex_costs / np.max(ex_costs) >= slack].max()
        ]
        ex_thr_ns = np.sort(_ex)[
            np.arange(len(_ex))[ex_costs / np.max(ex_costs) >= 1].max()
        ]
        return ex_thr, ex_thr_ns, ex_costs

    # Set K thresholds
    K_thrs = np.array(
        [
            getelbow(K, True),
            getelbow2(K, True),
            np.average(K, weights=1.0 / varex),
            fmin,
            fmid,
        ]
    )
    K_thr_safe, K_thr_agg = np.percentile(K_thrs, [25, 50], method="linear")

    # Determine dataset conditions
    high_grappa = False
    high_rho = False
    if (K.max() < fmax and len(rs) / nc > 0.75) or np.sum(
        e.zu_ex[: int(len(rs) * 0.15)]
    ) > 5:
        high_grappa = True
        # This if-condition shows the dataset to have low max K, yet K>R almost everywhere
        #  suggesting a dataset affectd by heavy GRAPPA/T2* weighted artifact
        #  so we raise the Kappa thresholds
        K_thr_safe, K_thr_agg = np.percentile(K_thrs, [50, 75], method="linear")

    # High Rho dataset
    if getelbow2(R, True) > 1.5 * fmin or np.mean(R) > 1.5 * fmin:
        high_rho = True
    R_thr = np.min([np.mean([getelbow(R, True), getelbow2(R, True)]), 1.5 * fmin])

    # Set optimization set
    if high_grappa or high_rho:
        # Determine ex thr from enriched set
        rs_opt = N[score > np.max([getelbow2(score, True), np.log10(score.max())])]
        rs_opt = np.setdiff1d(rs_opt, N[zu_ex > 5])
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
    acc_a = rs[acc_sel_a == 3]
    eject_a_ss0 = acc_a[
        rankvec(ss0[acc_a]) - rankvec(K[acc_a]) >= np.max([6, len(acc_a) / 2])
    ]
    eject_a_R = acc_a[
        rankvec(R[acc_a]) - rankvec(K[acc_a]) >= np.max([6, len(acc_a) / 2])
    ]
    eject_a_er = acc_a[
        rankvec(er[acc_a]) - rankvec(K[acc_a]) >= np.max([6, len(acc_a) / 2])
    ]
    eject_a = np.union1d(eject_a_R, eject_a_ss0)
    eject_a = np.union1d(eject_a, eject_a_er)
    acc_a = np.setdiff1d(acc_a, eject_a)
    acc_a = np.setdiff1d(acc_a, acc_a[e.lz[acc_a] == 0])

    # Nonparametric thresholds
    sr2_ex = sr2 / ex
    sr2_ss0 = sr2 - ss0
    sr2_ex_thr = sr2_ex[acc_a].min() * 0.2
    sr2_ss0_thr = sr2_ss0[acc_a].min() * 0.2
    sr2_thr = np.percentile(sr2[acc_a], 20) * 0.2
    lz_thr = e.lz[acc_a].min()

    acc_sel_b = andb(
        [
            sr2[rs] > sr2_thr,
            sr2_ex[rs] > sr2_ex_thr,
            sr2_ss0[rs] > sr2_ss0_thr,
        ]
    )
    acc_b = rs[andb([acc_sel_b == 3, e.lz[rs] >= lz_thr]) > 0]
    acc_b = np.intersect1d(acc_b, rs[K[rs] > K_thr_safe])

    # Clean up based on R scores in acc_b while protecting robust R2* task effects
    z_sr2_ss0_a = ((sr2 - ss0)[N] - (sr2 - ss0)[acc_a].mean()) / (sr2 - ss0)[
        acc_a
    ].std()
    # eject_R = acc_b[z_sr2_ss0_a[acc_b] > 4]
    z_R_a = (R[N] - R[acc_a].mean()) / R[acc_a].std()
    z_ss0_a = (ss0[N] - ss0[acc_a].mean()) / (ss0)[acc_a].std()
    eject_R_hiK = acc_b[
        andb(
            [
                (z_R_a - z_sr2_ss0_a)[acc_b] > 0,
                z_R_a[acc_b] > 0,
                z_ss0_a[acc_b] > 3,
                N[acc_b] < np.min([getelbow(K), 15]),
            ]
        )
        == 4
    ]
    eject_R_lowK = acc_b[
        andb(
            [
                (z_R_a - z_sr2_ss0_a)[acc_b] > 0,
                z_R_a[acc_b] > 0,
                # z_ss0_a[acc_b] > 0,
                N[acc_b] >= np.min([getelbow(K), 15]),
            ]
        )
        == 3
    ]
    eject_R_highS0 = acc_b[z_ss0_a[acc_b] > 15]
    eject_R = np.union1d(np.union1d(eject_R_hiK, eject_R_lowK), eject_R_highS0)

    # Remove based on high relative ex compared to acc_a
    # import pudb
    #
    # pudb.set_trace()

    z_ex_a = (ex[N] - ex[acc_a].mean()) / ex[acc_a].std()
    z_sr2_ex_a = (sr2_ex[N] - sr2_ex[acc_a].mean()) / sr2_ex[acc_a].std()
    # zu_lz = zu(e.lz)
    eject_ex_hiK = acc_b[
        andb(
            [
                (z_ex_a - z_sr2_ex_a)[acc_b] > 6.5,
                z_ex_a[acc_b] > np.max([np.percentile(z_ex_a[K < K_thr_agg], 20), 4]),
                N[acc_b] <= np.min([getelbow(K), 15]),
            ]
        )
        == 3
    ]
    eject_ex_lowK = acc_b[
        andb(
            [
                (z_ex_a - z_sr2_ex_a)[acc_b] > 6.5,
                z_ex_a[acc_b] > np.max([np.percentile(z_ex_a[K < K_thr_agg], 20), 4]),
                # z_ex_a[acc_b] > 4,
                N[acc_b] > np.min([getelbow(K), 15]),
                e.lz[acc_b] < np.percentile(e.lz[acc_a], 20),
            ]
        )
        == 4
    ]
    eject_ex = np.union1d(eject_ex_hiK, eject_ex_lowK)

    # Remove large included artifacts after tail inclusion
    eject_score_thr = np.min(
        [100, np.mean([getelbow(score, True), np.percentile(score[acc_a], 33)])]
    )
    high_varex_thr = np.percentile(varex[acc_a], 80)

    eject_quality = acc_b[
        andb(
            [
                score[acc_b] < eject_score_thr,
                R[acc_b] > R_thr,
                varex[acc_b] > high_varex_thr,
            ]
        )
        > 1
    ]
    eject_score_low = acc_b[score[acc_b] < np.log10(eject_score_thr)]
    eject_empty = acc_b[e.lz[acc_b] == 0]

    # High field artifact
    # e.emx = emx = np.max(em, axis=1)
    # e.emr = emr = (np.max(em, axis=1) / np.min(em, axis=1))
    eject_field_em = np.intersect1d(
        rs[z_(e.emx[rs]) > 9],
        rs[z_(e.emr[rs]) > 9],
    )
    eject_field_rank = acc_b[rankvec(e.emx[acc_b]) - rankvec(K[acc_b]) > len(acc_a)]
    eject_field = np.union1d(eject_field_em, eject_field_rank)
    # eject_field_agg = np.intersect1d(eject_field_em, acc_b[K[acc_b] < K_thrs.max()])
    eject_field_t2s = rs[e.zu_cT2s[rs] > 6.5]
    # eject_field_low_t2s = rs[
    #     andb([e.zu_cT2s[rs] > 6, score[rs] < eject_score_thr]) == 2
    # ]
    # eject_field_t2s = np.union1d(eject_field_low_t2s, eject_field_t2s)

    # Combine safe ejects
    eject_R = np.setdiff1d(eject_R, acc_a)  # R not allowed to go after acc_a
    eject = np.union1d(
        eject_field, eject_R
    )  # field will be protected by quality requirements
    eject = np.intersect1d(eject, eject_quality)
    eject = np.union1d(eject, eject_field_t2s)
    eject = np.union1d(eject, eject_empty)

    # Combine all ejects
    eject_ex = np.setdiff1d(eject_ex, acc_a)
    eject = np.union1d(eject, eject_ex)
    # eject = np.union1d(eject, eject_field_agg)
    eject = np.union1d(eject, eject_score_low)

    # Final acceptance
    acc = np.setdiff1d(acc_b, eject)

    min_score = score[acc].min()
    print("Minimum score is", min_score)

    rs = np.setdiff1d(rs, acc)
    # eject_ex_midk = np.setdiff1d(eject_ex, acc)
    # midk = np.union1d(rs[zu_ex[rs] > 5], eject_ex_midk)
    midk = rs[zu_ex[rs] > 5]

    ign = np.setdiff1d(rs, midk)
    ign = np.union1d(ign, eject)

    e.C = np.ones(nc) * 9
    e.C[acc] = 1
    e.C[rej] = 0
    e.C[ign] = 3
    e.C[eject] = 2
    e.C[midk] = -1

    return e
