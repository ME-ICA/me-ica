import gzip
import os
import pickle

import numpy as np

from .fastica import FastICA
from .fitmodels import fitmodels_direct, fitmodels_pca
from .utils.selection import andb, getelbow, getfbounds
from .utils.volume import eimask, fmask, makemask, unmask


def min_dFdV(ctpca, krsel=1):
    """
    In this function we maximimze the L1 norm of the
    derivative of variance versus TE-dependence, i.e.
    F s.t. max(|dV/dF|_1), so essentially an
    "L1 derivative", which is done in signal processesing
    for tasks such as edge detection.
    """
    assert krsel in [1, 2]
    ctpca_v = ctpca[ctpca[:, 3].argsort()[::-1], :]
    kr = np.log(ctpca[:, krsel])
    nstep = kr.max() / 0.01
    ls_dF = np.linspace(kr.min(), kr.max(), int(nstep))
    assert ls_dF.max() == kr.max()
    ndFdV0 = []
    for _i, _F in enumerate(ls_dF):
        # kapaps
        plsel = (ctpca_v[:, krsel] > np.exp(_F)).astype(int)
        sel_grad = np.gradient(plsel)
        ndFdV0.append(np.count_nonzero(sel_grad))
    F_thr = np.exp(np.array(ls_dF)[np.argmax(ndFdV0)])
    return F_thr


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

    # if wvpca:
    #       import ipdb
    #       ipdb.set_trace()
    #       print "++Transforming time series from time domain to wavelet domain."
    #       dz,cAl = dwtmat(dz)

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

        # ipdb.set_trace()

        eimum = np.atleast_2d(eim)
        eimum = np.transpose(eimum, np.argsort(np.atleast_2d(eim).shape)[::-1])
        eimum = np.array(np.squeeze(unmask(eimum.prod(1), assets.mask)), dtype=bool)
        vTmix = v.T
        vTmixN = ((vTmix.T - vTmix.T.mean(0)) / vTmix.T.std(0)).T
        # ctb,KRd,betasv,v_T = fitmodels2(catd,v.T,eimum,t2s,tes,mmixN=vTmixN)

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

    # import pudb; pudb.set_trace()

    # Select a high dimensionli subspace
    fmin, fmid, fmax = getfbounds(ne)
    kappa_elbow = getelbow(ctb[:, 1], True)
    kappa_thr_det = min_dFdV(ctb)
    kmin = ctb[:, 1].min()
    kappa_thr = np.average(
        sorted([kappa_thr_det, kmin, kappa_elbow]),
        weights=(assets.kdaw, 1, 1),  # type: ignore
    )
    pcsel_kappa = ctb[:, 1] > kappa_thr

    rho_thr_det = min_dFdV(ctb, 2)
    rho_elbow = getelbow(ctb[:, 2], True)
    rmin = ctb[:, 2].min()
    rho_thr = np.average(
        sorted([rho_thr_det, rmin, rho_elbow]),
        weights=(assets.rdaw, 1, 1),  # type: ignore
    )
    pcsel_rho = ctb[:, 2] > rho_thr

    # Make sure we keep top-most variance PCs
    # import pudb; pudb.set_trace()
    pcsel_var = ctb[:, -1] > eigelb

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
    return mmix
