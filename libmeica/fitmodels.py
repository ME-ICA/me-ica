import os
from datetime import datetime
from typing import List, Optional, cast

import numpy as np
from scipy.stats import spearmanr, zscore

from libmeica.t2smap import optcom
from libmeica.utils.lstsq import get_coeffs

from .utils.filter import spatclust
from .utils.memory import create_memmap, flush_after
from .utils.selection import getfbounds, mmix_hash, save_variables_to_file, z_
from .utils.volume import fmask, niwrite, unmask

F_MAX = 500
Z_MAX = 8


def computefeats2(data, mmix, mask, normalize=True, assets=None):
    # Write feature versions of components
    data = data[mask]
    data_vn = (data - data.mean(axis=-1)[:, np.newaxis]) / data.std(axis=-1)[
        :, np.newaxis
    ]
    if assets is None:
        v_nosel = None
    else:
        v_nosel = assets.v_nosel
    data_R = get_coeffs(unmask(data_vn, mask), mask, mmix, add_const=v_nosel)[mask]
    data_R[data_R < -0.999] = -0.999
    data_R[data_R > 0.999] = 0.999
    data_Z = np.arctanh(data_R)
    if len(data_Z.shape) == 1:
        data_Z = np.atleast_2d(data_Z).T
    if normalize:
        # data_Z2 = ((data_Z.T-data_Z.mean(0)[:,np.newaxis])/data_Z.std(0)[:,np.newaxis]).T
        data_Z = (
            ((data_Z.T - data_Z.mean(0)[:, np.newaxis]) / data_Z.std(0)[:, np.newaxis])
            + (data_Z.mean(0) / data_Z.std(0))[:, np.newaxis]
        ).T
    return data_Z


def fitmodels_pca(
    catd,
    mmix,
    mask,
    t2s,
    tes,
    *,
    reindex=False,
    mmixN=None,
    assets=None,
):
    """
    Usage:

    fitmodels_pca(fout)

    Input:
    fout is flag for output of per-component TE-dependence maps
    t2s is a (nx,ny,nz) ndarray
    tes is a 1d array
    """

    # Compute opt. com. raw data
    tsoc = assets.OCcatd[mask]
    # tsoc = np.array(optcom(catd, t2s, tes, mask), dtype=float)[mask]
    tsoc_mean = tsoc.mean(axis=-1)
    tsoc_dm = tsoc - tsoc_mean[:, np.newaxis]

    # Compute un-normalized weight dataset (features)
    if mmixN is None:
        mmixN = mmix
    # WTS = computefeats2(unmask(unmask(tsoc,mask)[t2s!=0],t2s!=0),mmixN,t2s!=0,normalize=False)
    WTS = computefeats2(unmask(tsoc, mask), mmixN, mask, normalize=False)

    # Compute PSC dataset - shouldn't have to refit data
    # n.b. whole mixing matrix received here, so no add_const
    tsoc_B = get_coeffs(unmask(tsoc_dm, mask), mask, mmix)[mask]

    totvar = (tsoc_B**2).sum()
    totvar_norm = (WTS**2).sum()

    # Compute Betas and means over TEs for TE-dependence analysis
    Ne = tes.shape[0]
    nx, ny, nz, _, nt = catd.shape
    nc = mmix.shape[1]
    betas = create_memmap(
        "_pca_betas",
        shape=(
            nx,
            ny,
            nz,
            Ne,
            nc,
        ),
    )
    betas[:] = np.zeros([nx, ny, nz, Ne, nc], dtype=np.float32)
    for _e in range(Ne):
        catd_e = catd[:, :, :, _e, :].copy()
        get_coeffs(catd_e, mask, mmix, betas5d=betas, Ne5d=_e)
    nx, ny, nz, Ne, nc = betas.shape
    NmD = (t2s != 0).sum()
    mu = catd.mean(axis=-1)
    tes = np.reshape(tes, (Ne, 1))

    # Mask arrays
    mumask = fmask(mu, t2s != 0)
    betamask = fmask(betas, t2s != 0)

    # Setup Xmats
    # Model 1
    X1 = mumask.transpose()

    # Model 2
    X2 = -1 * np.tile(tes, (1, NmD)) * mumask.transpose()  # type: ignore

    # Tables for component selection
    Kappas = np.zeros([nc])
    Rhos = np.zeros([nc])
    varex = np.zeros([nc])
    varex_norm = np.zeros([nc])

    print("Evaluating TE-Dependence of PCA", end="", flush=True)

    for i in range(nc):
        # size of B is (nc, nx*ny*nz)
        # with flush_after(*mmap_list):
        B = np.atleast_3d(betamask)[:, :, i].transpose()
        alpha = (np.abs(B) ** 2).sum(axis=0)
        varex[i] = (tsoc_B[:, i] ** 2).sum() / totvar * 100.0
        varex_norm[i] = (
            (unmask(WTS, mask)[t2s != 0][:, i] ** 2).sum() / totvar_norm * 100.0
        )

        # S0 Model
        coeffs_S0 = (B * X1).sum(axis=0) / (X1**2).sum(axis=0)
        SSE_S0 = (B - X1 * np.tile(coeffs_S0, (Ne, 1))) ** 2  # type: ignore
        SSE_S0 = SSE_S0.sum(axis=0)
        F_S0 = (alpha - SSE_S0) * 2 / (SSE_S0)

        # R2 Model
        coeffs_R2 = (B * X2).sum(axis=0) / (X2**2).sum(axis=0)
        SSE_R2 = (B - X2 * np.tile(coeffs_R2, (Ne, 1))) ** 2
        SSE_R2 = SSE_R2.sum(axis=0)
        F_R2 = (alpha - SSE_R2) * 2 / (SSE_R2)

        # Compute weights as Z-values
        wtsZ = (WTS[:, i] - WTS[:, i].mean()) / WTS[:, i].std()
        wtsZ[np.abs(wtsZ) > Z_MAX] = (Z_MAX * (np.abs(wtsZ) / wtsZ))[
            np.abs(wtsZ) > Z_MAX
        ]

        # Compute Kappa and Rho
        F_S0[F_S0 > F_MAX] = F_MAX
        F_R2[F_R2 > F_MAX] = F_MAX
        F_S0[np.isnan(F_S0)] = 0
        F_R2[np.isnan(F_R2)] = 0
        F_S0[np.isinf(F_S0)] = 0
        F_R2[np.isinf(F_R2)] = 0
        Kappas[i] = np.average(
            F_R2, weights=np.abs(np.squeeze(unmask(wtsZ, mask)[t2s != 0] ** 2.0))
        )
        Rhos[i] = np.average(
            F_S0, weights=np.abs(np.squeeze(unmask(wtsZ, mask)[t2s != 0] ** 2.0))
        )
        print(".", end="", flush=True)

    comptab_pre = np.vstack([np.arange(nc), Kappas, Rhos, varex, varex_norm]).T
    if reindex:
        # Re-index all components in Kappa order
        comptab = comptab_pre[comptab_pre[:, 1].argsort()[::-1], :]
        nnc = np.array(comptab[:, 0], dtype=int)
        mmix_new = mmix[:, nnc]
        comptab[:, 0] = np.arange(comptab.shape[0])
    else:
        comptab = comptab_pre
        mmix_new = mmix
    return comptab, betas, mmix_new


def fitmodels_direct(
    catd,
    mmix,
    mask,
    t2s,
    tes,
    *,
    fout=None,
    reindex=False,
    mmixN=None,
    full_sel=True,
    debugout=False,
    assets=None,
):
    """
    Usage:

    fitmodels_direct(fout)

    Input:
    fout is flag for output of per-component TE-dependence maps
    t2s is a (nx,ny,nz) ndarray
    tes is a 1d array
    """

    if assets is not None:
        aff = assets.aff
        head = assets.head
        ne = assets.ne
        args = assets.args
        v_nosel = assets.v_nosel
    else:
        v_nosel = None

    # Compute opt. com. raw data
    tsoc = np.array(optcom(catd, t2s, tes, mask), dtype=float)[mask]
    tsoc_mean = tsoc.mean(axis=-1)
    tsoc_dm = tsoc - tsoc_mean[:, np.newaxis]

    # Compute un-normalized weight dataset (features)
    if mmixN is None:
        mmixN = mmix
    # WTS = computefeats2(unmask(unmask(tsoc,mask)[t2s!=0],t2s!=0),mmixN,t2s!=0,normalize=False)
    WTS = computefeats2(unmask(tsoc, mask), mmixN, mask, normalize=False, assets=assets)

    # Compute PSC dataset - shouldn't have to refit data
    tsoc_B = get_coeffs(unmask(tsoc_dm, mask), mask, mmix, add_const=v_nosel)[mask]
    tsoc_Babs = np.abs(tsoc_B)
    PSC = tsoc_B / tsoc.mean(axis=-1)[:, np.newaxis] * 100

    # Compute skews to determine signs based on unnormalized weights, correct mmix & WTS signs based on spatial distribution tails
    from scipy.stats import skew

    signs = skew(WTS, axis=0)
    signs /= np.abs(signs)
    mmix = mmix.copy()
    mmix *= signs
    WTS *= signs
    PSC *= signs
    totvar = (tsoc_B**2).sum()
    totvar_norm = (WTS**2).sum()

    # Compute Betas and means over TEs for TE-dependence analysis
    Ne = tes.shape[0]

    nx, ny, nz, _, nt = catd.shape
    nc = mmix.shape[1]
    print("nc is", nc, mmix.shape)
    betas = np.zeros([nx, ny, nz, Ne, nc], dtype=np.float32)
    for _e in range(Ne):
        get_coeffs(
            catd[:, :, :, _e, :],
            mask,
            mmix,
            betas5d=betas,
            Ne5d=_e,
            add_const=v_nosel,
        )
    # betas = np.array(betas).transpose([1,2,3,0,4])
    # betas = cat2echos(
    #     get_coeffs(uncat2echos(catd, Ne), np.tile(mask, (1, 1, Ne)), mmix), Ne
    # )
    nx, ny, nz, Ne, nc = betas.shape
    Nm = mask.sum()
    NmD = (t2s != 0).sum()
    mu = catd.mean(axis=-1)
    tes = np.reshape(tes, (Ne, 1))
    fmin, fmid, fmax = getfbounds(ne)  # type: ignore
    csize = np.max([int(Nm * 0.0005) + 5, 20])

    # Mask arrays
    mumask = fmask(mu, t2s != 0)
    t2smask = fmask(t2s, mask)
    t2smask = fmask(t2s, t2s != 0)
    betamask = fmask(betas, t2s != 0)

    # PSC by TE
    psc_by_te_mask = (betamask / mumask[:, :, np.newaxis]) * 100

    if debugout:
        fout = aff  # type: ignore

    # Setup Xmats
    # Model 1
    X1 = mumask.transpose()

    # Model 2
    X2 = -1 * np.tile(tes, (1, NmD)) * mumask.transpose()  # type: ignore

    # Tables for component selection
    Kappas = np.zeros([nc])
    Rhos = np.zeros([nc])
    varex = np.zeros([nc])
    corrT1mi = np.zeros([nc])
    corrT2s = np.zeros([nc])
    varex_norm = np.zeros([nc])
    F_R2_lclsz = np.zeros([nc])
    F_S0_lclsz = np.zeros([nc])
    Z_lclsz = np.zeros([nc])

    # Assuming Nm, NmD, nc are already defined
    Z_maps = create_memmap("Z_maps", (Nm, nc))
    F_R2_maps = create_memmap("F_R2_maps", (NmD, nc))
    F_S0_maps = create_memmap("F_S0_maps", (NmD, nc))
    Z_clmaps = create_memmap("Z_clmaps", (Nm, nc))
    F_R2_clmaps = create_memmap("F_R2_clmaps", (NmD, nc))
    F_S0_clmaps = create_memmap("F_S0_clmaps", (NmD, nc))
    # Parametric outputs
    coeff_R2_maps = create_memmap("coeff_R2_maps", (NmD, nc))
    dS0_maps = create_memmap("dS0_maps", (NmD, nc))
    dT2_maps = create_memmap("dT2_maps", (NmD, nc))
    mmap_list = [
        Z_maps,
        F_R2_maps,
        F_S0_maps,
        Z_clmaps,
        F_R2_clmaps,
        F_S0_clmaps,
        coeff_R2_maps,
        dS0_maps,
        dT2_maps,
    ]

    # Z_maps = np.zeros([Nm, nc], dtype=np.float32)
    # F_R2_maps = np.zeros([NmD, nc], dtype=np.float32)
    # F_S0_maps = np.zeros([NmD, nc], dtype=np.float32)
    # Z_clmaps = np.zeros([Nm, nc], dtype=np.float32)
    # F_R2_clmaps = np.zeros([NmD, nc], dtype=np.float32)
    # F_S0_clmaps = np.zeros([NmD, nc], dtype=np.float32)
    # # Br_clmaps_R2 = np.zeros([Nm, nc])
    # # Br_clmaps_S0 = np.zeros([Nm, nc])

    print("Evaluating TE-Dependence of ICA components", end="", flush=True)

    for i in range(nc):
        # size of B is (nc, nx*ny*nz)
        # with flush_after(*mmap_list):
        B = np.atleast_3d(betamask)[:, :, i].transpose()
        alpha = (np.abs(B) ** 2).sum(axis=0)
        varex[i] = (tsoc_B[:, i] ** 2).sum() / totvar * 100.0
        varex_norm[i] = (
            (unmask(WTS, mask)[t2s != 0][:, i] ** 2).sum() / totvar_norm * 100.0
        )

        # S0 Model
        coeffs_S0 = (B * X1).sum(axis=0) / (X1**2).sum(axis=0)
        SSE_S0 = (B - X1 * np.tile(coeffs_S0, (Ne, 1))) ** 2  # type: ignore
        SSE_S0 = SSE_S0.sum(axis=0)
        F_S0 = (alpha - SSE_S0) * 2 / (SSE_S0)
        F_S0_maps[:, i] = F_S0
        dS0_maps[:, i] = coeffs_S0

        # R2 Model
        # import ipdb
        # ipdb.set_trace()
        coeffs_R2 = (B * X2).sum(axis=0) / (X2**2).sum(axis=0)
        SSE_R2 = (B - X2 * np.tile(coeffs_R2, (Ne, 1))) ** 2
        SSE_R2 = SSE_R2.sum(axis=0)
        F_R2 = (alpha - SSE_R2) * 2 / (SSE_R2)
        F_R2_maps[:, i] = F_R2
        coeff_R2_maps[:, i] = coeffs_R2  # coeffs_R2 are delta-R2*

        # Derive dT2* from dR2*   # testing
        # import pdb
        # pdb.set_trace()
        # dT2_maps[:, i] = 1.0 / ((1 / t2smask) + coeffs_R2) - t2smask
        # pdT2_maps[:, i] = dT2_maps[:, i] / t2smask * 100

        # Compute weights as Z-values
        wtsZ = (WTS[:, i] - WTS[:, i].mean()) / WTS[:, i].std()
        wtsZ[np.abs(wtsZ) > Z_MAX] = (Z_MAX * (np.abs(wtsZ) / wtsZ))[
            np.abs(wtsZ) > Z_MAX
        ]
        Z_maps[:, i] = wtsZ

        _F_R2_clmap = np.squeeze(unmask(F_R2_maps[:, i], t2s != 0))
        _F_S0_clmap = np.squeeze(unmask(F_S0_maps[:, i], t2s != 0))
        _Z_clmap = np.squeeze(unmask(Z_maps[:, i], mask))

        _clmap, _clszs = spatclust(_F_R2_clmap, thr=fmin, csize=csize, getsizes=True)
        F_R2_clmaps[:, i] = fmask(_clmap != 0, t2s != 0)
        if 0 in _clszs.keys() and len(_clszs.items()) > 1:
            del _clszs[0]
            F_R2_lclsz[i] = np.max(list(_clszs.values()))

        _clmap, _clszs = spatclust(_F_S0_clmap, thr=fmin, csize=csize, getsizes=True)
        F_S0_clmaps[:, i] = fmask(_clmap != 0, t2s != 0)
        if 0 in _clszs.keys() and len(_clszs.items()) > 1:
            del _clszs[0]
            F_S0_lclsz[i] = np.max(list(_clszs.values()))

        _clmap, _clszs = spatclust(_Z_clmap, thr=1.95, csize=csize, getsizes=True)
        Z_clmaps[:, i] = fmask(_clmap != 0, mask)
        if 0 in _clszs.keys() and len(_clszs.items()) > 1:
            del _clszs[0]
            Z_lclsz[i] = np.max(list(_clszs.values()))

        # Compute Kappa and Rho
        F_S0[F_S0 > F_MAX] = F_MAX
        F_R2[F_R2 > F_MAX] = F_MAX
        Kappas[i] = np.average(
            F_R2, weights=np.abs(np.squeeze(unmask(wtsZ, mask)[t2s != 0] ** 2.0))
        )
        Rhos[i] = np.average(
            F_S0, weights=np.abs(np.squeeze(unmask(wtsZ, mask)[t2s != 0] ** 2.0))
        )

        # Calculate correlation to MIR image and T2* map
        t1mimask = fmask(assets.t1mi, mask)
        corrT1mi_sp = spearmanr(t1mimask, WTS[:, i])
        corrT1mi[i] = corrT1mi_sp.correlation
        t2smask_z = zscore(t2smask)

        wtsZ_t2s = wtsZ.copy()
        wtsZ_t2s[wtsZ_t2s < -4] = -4
        wtsZ_t2s[wtsZ_t2s > 4] = 4
        t2s_z = np.average(t2smask_z, weights=wtsZ_t2s**2) * 100
        corrT2s[i] = t2s_z

        print(".", end="", flush=True)

    # Tabulate component values
    comptab_pre = np.vstack([np.arange(nc), Kappas, Rhos, varex, varex_norm]).T
    mmix_new = mmix
    if reindex:
        # Re-index all components in Kappa order
        comptab = comptab_pre[comptab_pre[:, 1].argsort()[::-1], :]
        Kappas = comptab[:, 1]
        Rhos = comptab[:, 2]
        varex = comptab[:, 3]
        varex_norm = comptab[:, 4]
        nnc = np.array(comptab[:, 0], dtype=int)
        mmix_new = mmix[:, nnc]
        F_S0_maps = F_S0_maps[:, nnc]
        F_R2_maps = F_R2_maps[:, nnc]
        Z_maps = Z_maps[:, nnc]
        WTS = WTS[:, nnc]
        PSC = PSC[:, nnc]
        tsoc_B = tsoc_B[:, nnc]
        tsoc_Babs = tsoc_Babs[:, nnc]
        coeff_R2_maps = coeff_R2_maps[:, nnc]
        dS0_maps = dS0_maps[:, nnc]
        dT2_maps = dT2_maps[:, nnc]
        psc_by_te_mask = psc_by_te_mask[..., nnc]
        F_R2_clmaps = F_R2_clmaps[..., nnc]
        F_S0_clmaps = F_S0_clmaps[..., nnc]
        F_R2_lclsz = F_R2_lclsz[..., nnc]
        F_S0_lclsz = F_S0_lclsz[..., nnc]
        Z_clmaps = Z_clmaps[..., nnc]
        Z_lclsz = Z_lclsz[..., nnc]
        comptab[:, 0] = np.arange(comptab.shape[0])
        corrT1mi = corrT1mi[nnc]
        corrT2s = corrT2s[nnc]
    else:
        comptab = comptab_pre
        mmix_new = mmix
    flush_after(*mmap_list)

    fit_timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    fit_hash = mmix_hash(timestring=fit_timestamp, mmix=mmix_new)
    mmix_id = f"{fit_hash}_{fit_timestamp}"
    if assets is not None:
        assets.mmix_id = mmix_id

    # Write out prelimary post-ICA comptab
    if full_sel:
        np.savetxt("comptab_ica_pre.txt", comptab, fmt="%.02f", delimiter="\t")

    # Full selection including clustering criteria
    seldict = {}
    selvars = [
        "nc",
        "Kappas",
        "Rhos",
        "WTS",
        "varex",
        "Z_maps",
        "F_R2_maps",
        "F_S0_maps",
        "Z_clmaps",
        "F_R2_clmaps",
        "F_S0_clmaps",
        "F_R2_lclsz",
        "F_S0_lclsz",
        "Z_lclsz",
        "tsoc_B",
        "PSC",
        "mmix_new",
        "Ne",
        "mask",
        "t2s",
        "head",
        "dS0_maps",  # for parametric outputs
        "dT2_maps",  # and below
        "psc_by_te_mask",
        "coeff_R2_maps",
        "corrT1mi",
        "corrT2s",
        "aff",
        "mmix_id",
    ]
    for vv in selvars:
        seldict[vv] = eval(vv)

    if debugout or ("DEBUGOUT" in list(args)):  # type: ignore
        debugoutvars = {
            "seldict": seldict,
            "Ne": Ne,
            "mask": mask,
            "t2s": t2s,
            "head": head,  # type: ignore
            "aff": aff,  # type: ignore
        }
        save_variables_to_file(debugoutvars, "selection_debug.pkl.gz")

    # if fout_all:
    #     seldict_to_fout(seldict, assets=assets)
    return seldict, comptab, betas, mmix_new


def seldict_to_fout(seldict, *, assets, component_list: Optional[List] = None):
    PSC = seldict["PSC"]
    F_R2_maps = seldict["F_R2_maps"]
    F_S0_maps = seldict["F_S0_maps"]
    Z_maps = seldict["Z_maps"]
    dS0_maps = seldict["dS0_maps"]
    dT2_maps = seldict["dT2_maps"]
    psc_by_te_mask = seldict["psc_by_te_mask"]
    coeff_R2_maps = seldict["coeff_R2_maps"]
    # nc = seldict["nc"]

    if component_list is None:
        component_list = cast(list, assets.acc)

    for i in component_list:
        # Save out files
        nx, ny, nz, nt, _ = assets.catd.shape
        out = np.zeros((nx, ny, nz, 9 + assets.Ne))
        ccname = f"cc{i:03d}.nii"

        out[:, :, :, 0] = np.squeeze(unmask(PSC[:, i], assets.mask))
        out[:, :, :, 1] = np.squeeze(unmask(F_R2_maps[:, i], assets.t2s != 0))
        out[:, :, :, 2] = np.squeeze(unmask(F_S0_maps[:, i], assets.t2s != 0))
        out[:, :, :, 3] = np.squeeze(unmask(Z_maps[:, i], assets.mask))
        out[:, :, :, 4] = np.squeeze(
            unmask(coeff_R2_maps[:, i], assets.t2s != 0)
        )  # delta-R2* -> subbrik 4
        out[:, :, :, 5] = np.squeeze(unmask(dS0_maps[:, i], assets.t2s != 0))
        out[:, :, :, 6] = np.squeeze(unmask(dT2_maps[:, i], assets.t2s != 0))

        for _psc_te in range(assets.Ne):
            # from pudb import set_trace; set_trace()
            out[:, :, :, 9 + _psc_te] = np.squeeze(
                unmask(psc_by_te_mask[:, _psc_te, i], assets.t2s != 0)
            )

        niwrite(out, assets.aff, ccname, header=assets.head)
        psc_labelstrings = " ".join(
            [f"-sublabel {i + 8 + 1} PSC_te{i + 1}" for i in range(assets.Ne)]
        )
        os.system(
            "3drefit -sublabel 0 PSC -sublabel 1 F_R2  -sublabel 2 F_SO -sublabel 3 Z_sn -sublabel 4 dR2s  -sublabel 5 dS0 %s %s 2>  /dev/null > /dev/null"
            % (psc_labelstrings, ccname)
        )
