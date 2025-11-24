import gzip
import os
import pickle

import numpy as np

from .fastica import FastICA
from .fitmodels import fitmodels_direct, fitmodels_pca
from .utils.selection import andb, getelbow, getelbow2, getfbounds
from .utils.volume import eimask, fmask, makemask, unmask
from .utils.filter import rankvec


def skew_positive(x, alpha=0.5, n=None):
    """
    Smoothly boost positive values of x.

    x     : scalar or array
    alpha : how much to increase + side (0 = no skew)
    n     : roughly max abs value in your data;
            if None, it's inferred from x
    """
    x = np.asarray(x, dtype=float)
    if n is None:
        n = np.max(np.abs(x)) if np.max(np.abs(x)) > 0 else 1.0

    k = 4.0 / n  # width of transition region
    sigma = 1.0 / (1.0 + np.exp(-k * x))
    return x * (1.0 + alpha * sigma)


def KR_V_thr(K, V, daw):
    assert daw > 0
    drank = rankvec(K)[::-1] - rankvec(V)[::-1]
    if daw > 1 and daw <= 10:
        thr = np.average(
            K,
            weights=(1.1) ** skew_positive(drank, alpha=1 + (daw / 10)),
        )
    elif daw <= 1:
        V_weight = 1 - daw
        K_weight = daw
        K_thr_V = np.average(K, weights=V)
        K_thr_rank = np.average(
            K,
            weights=(1.1) ** skew_positive(drank, alpha=1.1),
        )
        thr = np.average([K_thr_V, K_thr_rank], weights=[V_weight, K_weight])
    else:  # daw > 10
        min_weight = daw / 10
        K_weight = 1
        K_thr = np.average(
            K,
            weights=(1.1) ** skew_positive(drank, alpha=1 + (10 / 10)),
        )
        K_min = K.min()
        thr = np.average(
            [K_thr, K_min],
            weights=[K_weight, min_weight],
        )
    return thr


def KR_V_thr_test(K, V):
    daws = np.linspace(0.01, 100, 500)
    daw_res = []
    for daw in daws:
        daw_res.append([daw, np.sum(K > KR_V_thr(K, V, daw))])
    np.savetxt("daw_test.1D", daw_res, delimiter=" ")


def tedpca(
    ste=0,
    *,
    assets,
    mlepca=True,
):
    catd = assets.catd
    nx, ny, nz, ne, nt = catd.shape
    ste = np.array([int(ee) for ee in str(ste).split(",") if ee != ""])
    cAl = None
    if len(ste) == 1 and ste[0] == -1:
        print("-Computing PCA of optimally combined multi-echo data")
        OCmask = makemask(assets.OCcatd[:, :, :, np.newaxis, :])
        d = fmask(assets.OCcatd, OCmask)
        eim = eimask(d[:, np.newaxis, :])
        eim = eim[:, 0] == 1
        d = d[eim, :]
        # ipdb.set_trace()
    elif len(ste) == 1 and ste[0] == 0:
        print("-Computing PCA of spatially concatenated multi-echo data")
        ste = np.arange(ne)
        d = np.float64(fmask(catd, assets.mask))
        eim = eimask(d) == 1
        d = d[eim]  # type: ignore
    else:
        print("-Computing PCA of TE #%s" % ",".join([str(ee) for ee in ste]))
        d = np.float64(
            np.concatenate(
                [
                    fmask(catd[:, :, :, ee, :], assets.mask)[:, np.newaxis, :]
                    for ee in ste - 1
                ],
                axis=1,
            )
        )
        eim = eimask(d) == 1
        eim = np.squeeze(eim)
        d = np.squeeze(d[eim])  # type: ignore

    dz = ((d.T - d.T.mean(0)) / d.T.std(0)).T  # Variance normalize timeseries
    dz = (dz - dz.mean()) / dz.std()  # Variance normalize everything

    pcastate_fn = "pcastate.pklgz"

    if not os.path.exists(pcastate_fn):
        if mlepca:
            from sklearn.decomposition import PCA

            ppca = PCA(
                n_components=nt - 2,
                svd_solver="randomized",
                random_state=0,
                iterated_power=10,
            )
            ppca.fit(dz)
            v = ppca.components_
            s = ppca.explained_variance_
            u = np.dot(np.dot(dz, v.T), np.diag(1.0 / s))
        else:
            u, s, v = np.linalg.svd(dz, full_matrices=False)

        sp = s / s.sum()
        eigelb = sp[getelbow(sp)]

        spdif = np.abs(sp[1:] - sp[:-1])
        spdifh = spdif[spdif.shape[0] // 2 :]
        spdmin = spdif.min()
        spdthr = np.mean([spdifh.max(), spdmin])
        spmin = sp[
            (spdif.shape[0] // 2)
            + (np.arange(spdifh.shape[0])[spdifh >= spdthr][0])
            + 1
        ]
        spcum = []
        spcumv = 0
        for sss in sp:
            spcumv += sss
            spcum.append(spcumv)
        spcum = np.array(spcum)

        # Compute K and Rho for PCA comps

        eimum = np.atleast_2d(eim)
        eimum = np.transpose(eimum, np.argsort(np.atleast_2d(eim).shape)[::-1])
        eimum = np.array(np.squeeze(unmask(eimum.prod(1), assets.mask)), dtype=bool)
        vTmix = v.T
        vTmixN = ((vTmix.T - vTmix.T.mean(0)) / vTmix.T.std(0)).T

        ctb, betasv, v_T = fitmodels_pca(
            catd,
            v.T,
            eimum,
            assets.t2s,
            assets.tes,
            mmixN=vTmixN,
            assets=assets,
        )
        ctb = ctb[ctb[:, 0].argsort(), :]
        ctb = np.vstack([ctb.T[0:3], sp]).T

        # Save state
        print("Saving PCA")
        pcastate = {
            "u": u,
            "s": s,
            "v": v,
            "ctb": ctb,
            "eigelb": eigelb,
            "spmin": spmin,
            "spcum": spcum,
        }
        try:
            pcastate_f = gzip.open(pcastate_fn, "wb")
            pickle.dump(pcastate, pcastate_f)
            pcastate_f.close()
        except BaseException:
            print("Could not save PCA solution!")

    else:
        print("Loading PCA")
        pcastate_f = gzip.open(pcastate_fn, "rb")
        pcastate = pickle.load(pcastate_f)  # A SimpleNameSpace
        u = pcastate["u"]
        s = pcastate["s"]
        v = pcastate["v"]
        ctb = pcastate["ctb"]
        eigelb = pcastate["eigelb"]
        spmin = pcastate["spmin"]
        spcum = pcastate["spcum"]
        pcastate_f.close()
    np.savetxt("comp_table_pca.txt", ctb[ctb[:, 1].argsort(), :][::-1])
    np.savetxt("mepca_mix.1D", v[ctb[:, 1].argsort()[::-1], :].T)

    _K = ctb[:, 1]
    _R = ctb[:, 2]
    _V = ctb[:, 3]

    kappa_thr = KR_V_thr(_K, _V, assets.kdaw)
    pcsel_kappa = _K > kappa_thr

    KR_V_thr_test(_K, _V)

    rho_thr = KR_V_thr(_R, _V, assets.rdaw)
    pcsel_rho = _R > rho_thr

    pcsel_var = ctb[:, -1] > eigelb

    # Final pc selection
    pcsel = andb([pcsel_kappa, pcsel_rho, pcsel_var]) > 0

    dd = u.dot(np.diag(s * np.array(pcsel, dtype=int))).dot(v)

    nc = s[pcsel].shape[0]
    print(pcsel)
    print(
        "--Selected %i components. Minimum Kappa=%0.2f Rho=%0.2f"
        % (nc, kappa_thr, rho_thr)
    )

    dd = ((dd.T - dd.T.mean(0)) / dd.T.std(0)).T  # Variance normalize timeseries
    dd = (dd - dd.mean()) / dd.std()  # Variance normalize everything

    return nc, dd


def tedica(dd, nc, *, cost, assets):
    """
    Input is dimensionally reduced spatially concatenated multi-echo time series dataset from tedpca()
    Output is comptable, mmix, smaps from ICA, and betas from fitting catd to mmix
    """
    options = assets.options
    # Do ICA
    # from pudb import set_trace; set_trace()

    climit = float("%s" % options.conv)
    icanode = FastICA(
        n_components=nc,
        algorithm="parallel",
        fun="logcosh",
        tol=climit,
        max_iter=1000,
        whiten="unit-variance",
    )
    smaps = icanode.fit_transform(dd)
    mmix = np.array(icanode.mixing_)
    mmix = (mmix - mmix.mean(0)) / mmix.std(0)

    converge_success = icanode.n_iter_ <= climit

    return mmix, converge_success
