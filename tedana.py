#!/usr/bin/env python

import datetime
import os
import sys
import time

import nibabel as nib
import numpy as np

from libmeica.command_line import tedana_command_line, tedana_process_input
from libmeica.dr2s import signal_to_dr2s
from libmeica.fitmodels import (
    computefeats2,
    fitmodels_direct,
    get_coeffs,
    seldict_to_fout,
)
from libmeica.gscontrol import gscontrol_mmix, gscontrol_raw
from libmeica.mepca import tedica, tedpca
from libmeica.selcomps_encoding import SelcompsEncoding
from libmeica.t2smap import optcom
from libmeica.utils.memory import create_memmap
from libmeica.utils.t2s import t2sadmap
from libmeica.utils.volume import eimask, fmask, makeadmask, niwrite, unmask

try:
    from libmeica.vendor.output import vendor_outputs  # type: ignore

    VENDOR_MODE = True
except ImportError:
    VENDOR_MODE = False


__version__ = "v4.0.0"
welcome_block = """
# Multi-Echo ICA, Version %s
#
# Kundu, P., Brenowitz, N.D., Voon, V., Worbe, Y., Vertes, P.E., Inati, S.J., Saad, Z.S.,
# Bandettini, P.A. & Bullmore, E.T. Integrated strategy for improving functional
# connectivity mapping using multiecho fMRI. PNAS (2013).
#
# Kundu, P., Inati, S.J., Evans, J.W., Luh, W.M. & Bandettini, P.A. Differentiating
#   BOLD and non-BOLD signals in fMRI time series using multi-echo EPI. NeuroImage (2011).
# https://doi.org/10.1016/j.neuroimage.2011.12.028
#
# PROCEDURE 2 : Computes ME-PCA and ME-ICA
# -Computes T2* map
# -Computes PCA of concatenated ME data, then computes TE-dependence of PCs
# -Computes ICA of TE-dependence PCs
# -Identifies TE-dependent ICs, outputs high-Kappa (BOLD) component
#    and denoised time series
# -or- Computes TE-dependence of each component of a general linear model
#    specified by input (includes MELODIC FastICA mixing matrix)
""" % (__version__)


def write_split_ts(data, comptable, mmix, suffix="", *, glsig=None, assets=None):
    if assets is None:
        raise ValueError
    midk = assets.midk
    rej = assets.rej
    acc = assets.acc
    mask = assets.mask
    mdata = fmask(data, mask)
    v_nosel = assets.v_nosel
    # betas = fmask(
    #     get_coeffs(unmask((mdata.T - mdata.T.mean(0)).T, mask), mask, mmix), mask
    # )
    betas_tmp = mdata - mdata.mean(axis=1, keepdims=True)  # Demean: (V × T)
    betas_tmp = unmask(betas_tmp, mask)  # Unmask: (3D × T)
    betas_tmp = get_coeffs(
        betas_tmp, mask, mmix, add_const=v_nosel
    )  # Regression or fit: (3D × C)
    betas_tmp = fmask(betas_tmp, mask)  # Back to masked space: (V × C)
    betas = betas_tmp
    dmdata = mdata.T - mdata.T.mean(0)
    varexpl = (
        1 - ((dmdata.T - betas.dot(mmix.T)) ** 2.0).sum() / (dmdata**2.0).sum()
    ) * 100
    print("Variance explained: ", varexpl, "%")
    midkts = betas[:, midk].dot(mmix.T[midk, :])
    lowkts = betas[:, rej].dot(mmix.T[rej, :])

    betas_acc = betas[:, acc]

    if len(acc) != 0:
        _hik = betas_acc.dot(mmix.T[acc, :])
        if isinstance(glsig, np.ndarray):
            # import pudb; pudb.set_trace()
            _glsig = np.vstack([np.ones(glsig.shape[1]), glsig])
            # Does regression and projection in one step
            _hik_gscontrol = _hik - np.dot(
                np.linalg.lstsq(_glsig.T, _hik.T, rcond=None)[0].T, _glsig
            )
            niwrite(
                unmask(_hik_gscontrol, mask),
                assets.aff,
                f"hik_ts_{suffix}_T1c.nii",
                header=assets.head,
            )
        niwrite(
            unmask(_hik, mask),
            assets.aff,
            "_".join(["hik_ts", suffix]) + ".nii",
            header=assets.head,
        )
    # if len(midk) != 0:
    #     niwrite(
    #         unmask(midkts, mask),
    #         assets.aff,
    #         "_".join(["midk_ts", suffix]) + ".nii",
    #         header=assets.head,
    #     )
    # if len(rej) != 0:
    #     niwrite(
    #         unmask(lowkts, mask),
    #         assets.aff,
    #         "_".join(["lowk_ts", suffix]) + ".nii",
    #         header=assets.head,
    #     )
    niwrite(
        unmask(fmask(data, mask) - lowkts - midkts, mask),
        assets.aff,
        "_".join(["dn_ts", suffix]) + ".nii",
        header=assets.head,
    )
    return varexpl


def split_ts(
    *,
    data,
    assets,
):
    cbetas = get_coeffs(
        data - data.mean(-1)[:, :, :, np.newaxis],
        assets.mask,
        assets.mmix,
        add_const=assets.v_nosel,
    )
    betas = fmask(cbetas, assets.mask)
    if len(assets.acc) != 0:
        hikts = unmask(
            betas[:, assets.acc].dot(assets.mmix.T[assets.acc, :]), assets.mask
        )
    else:
        hikts = None
    return hikts, data - hikts


def writefeats(cbetas, comptable, mmix, *, assets, suffix=""):
    mask = assets.mask
    # Write signal changes (dS)
    niwrite(
        cbetas[:, :, :, :],
        assets.aff,
        "_".join(["betas", suffix]) + ".nii",
        header=assets.head,
    )
    niwrite(
        cbetas[:, :, :, assets.acc],
        assets.aff,
        "_".join(["betas_hik", suffix]) + ".nii",
        header=assets.head,
    )
    # Compute features (dS/S)
    if assets.options.e2d is None:
        e2d = np.floor(assets.ne / 2) + 1
    else:
        e2d = assets.options.e2d
    edm = fmask(assets.catd[:, :, :, e2d - 1, :], mask)
    edms = edm / edm.std(-1)[:, np.newaxis]
    edms[edm < 1] = 0
    hik, noise = split_ts(data=unmask(edms, mask), assets=assets)
    noise = noise - noise.mean(-1)[:, :, :, np.newaxis]

    zfac = (
        1.0 / (mmix.shape[0] - len(assets.acc) - 1) * (noise**2).sum(-1)
    )  # noise scaling
    niwrite(zfac, assets.aff, "zfac.nii", header=assets.head)

    cbetam = fmask(cbetas[:, :, :, assets.acc], mask)
    cbetam = (cbetam - cbetam.mean(0)) / cbetam.std(0)
    cbetam = cbetam / fmask(zfac, mask)[:, np.newaxis]
    cbetam[edm.mean(-1) < 1, :] = 0

    niwrite(
        unmask(cbetam, mask),
        assets.aff,
        "_".join(["feats", suffix]) + ".nii",
        header=assets.head,
    )


def writefeats2(data, mmix, mask, *, assets, suffix="", header=None):
    # Write feature versions of components
    if header is None:
        header = assets.head
    feats = computefeats2(data, mmix, mask)
    niwrite(
        unmask(feats, mask),
        assets.aff,
        "_".join(["feats", suffix]) + ".nii",
        header=assets.head,
    )


def writect(comptable, *, assets, ctname="", varexpl="-1", classarr=[]):
    if len(classarr) != 0:
        acc, rej, midk, empty = classarr
    else:
        acc = assets.acc
        midk = assets.midk
        rej = assets.rej
        empty = assets.empty
    nc = comptable.shape[0]
    ts = time.time()
    st = datetime.datetime.fromtimestamp(ts).strftime("%Y-%m-%d %H:%M:%S")
    sortab = comptable[comptable[:, 1].argsort()[::-1], :]
    if ctname == "":
        ctname = "comp_table.txt"
    open("accepted.txt", "w").write(",".join([str(int(cc)) for cc in acc]))
    open("rejected.txt", "w").write(",".join([str(int(cc)) for cc in rej]))
    open("midk_rejected.txt", "w").write(",".join([str(int(cc)) for cc in midk]))
    with open(ctname, "w") as f:
        f.write(
            "#\n#ME-ICA Component statistics table for: %s \n#Run on %s \n#\n"
            % (os.path.abspath(os.path.curdir), st)
        )
        f.write("#Dataset variance explained by ICA (VEx): %.02f \n" % (varexpl))
        f.write("#Total components generated by decomposition (TCo): %i \n" % (nc))
        f.write(
            "#No. accepted BOLD-like components, i.e. effective degrees of freedom for correlation (lower bound; DFe): %i\n"
            % (len(acc))
        )
        f.write(
            "#Total number of rejected components (RJn): %i\n" % (len(midk) + len(rej))
        )
        f.write(
            "#Nominal degress of freedom in denoised time series (..._medn.nii.gz; DFn): %i \n"
            % (assets.nt - len(midk) - len(rej))
        )
        f.write(
            "#ACC %s \t#Accepted BOLD-like components\n"
            % ",".join([str(int(cc)) for cc in acc])
        )
        f.write(
            "#REJ %s \t#Rejected non-BOLD components\n"
            % ",".join([str(int(cc)) for cc in rej])
        )
        f.write(
            "#MID %s \t#Rejected R2*-weighted artifacts\n"
            % ",".join([str(int(cc)) for cc in midk])
        )
        f.write(
            "#IGN %s \t#Ignored components (kept in denoised time series)\n"
            % ",".join([str(int(cc)) for cc in empty])
        )
        f.write("#VEx   TCo     DFe     RJn     DFn     \n")
        f.write(
            "##%.02f        %i      %i      %i      %i \n"
            % (
                varexpl,
                nc,
                len(acc),
                len(midk) + len(rej),
                assets.nt - len(midk) - len(rej),
            )
        )
        f.write("#      comp    Kappa   Rho     %%Var   %%Var(norm)     \n")
        for i in range(nc):
            f.write(
                "%d\t%f\t%f\t%.2f\t%.2f\n"
                % (sortab[i, 0], sortab[i, 1], sortab[i, 2], sortab[i, 3], sortab[i, 4])
            )


def writeresults(*, assets):
    comptable = assets.comptable
    mask = assets.mask
    print("++ Writing optimally combined time series")
    ts = assets.OCcatd
    niwrite(ts, assets.aff, "ts_OC.nii", header=assets.head)
    print("++ Writing Kappa-filtered optimally combined timeseries")
    varexpl = write_split_ts(ts, assets.comptable, assets.mmix, "OC", assets=assets)
    print("++ Writing component table")
    writect(comptable, ctname="comp_table.txt", varexpl=varexpl, assets=assets)
    print("++ Writing signal versions of components")
    ts_B = get_coeffs(ts, mask, assets.mmix, add_const=assets.v_nosel)
    if "DEBUG" in sys.argv and assets.mmix_id is not None:
        _suffix = "OC"
        mmix_suf = f"{_suffix}_{assets.mmix_id}"
        niwrite(
            ts_B[:, :, :, :],  # type: ignore
            assets.aff,
            "_".join(["betas", mmix_suf]) + ".nii",
            header=assets.head,
        )
    niwrite(
        ts_B[:, :, :, :],  # type: ignore
        assets.aff,
        "_".join(["betas", "OC"]) + ".nii",
        header=assets.head,
    )
    if len(assets.acc) != 0:
        niwrite(
            ts_B[:, :, :, assets.acc],  # type: ignore
            assets.aff,
            "_".join(["betas_hik", "OC"]) + ".nii",
            header=assets.head,
        )
        print("++ Writing optimally combined high-Kappa features")
        writefeats2(
            split_ts(data=ts, assets=assets)[0],
            assets.mmix[:, assets.acc],
            assets.mask,
            suffix="OC2",
            header=assets.head,
            assets=assets,
        )
    if assets.options.fout:
        # breakpoint()
        print("++ Writing TE-dependence SPMs of accepted components")
        seldict_to_fout(assets.seldict, assets=assets, component_list=assets.acc)


def writeresults_echoes(*, assets, glsig=None):
    for ii in range(assets.ne):
        print("++ Writing Kappa-filtered TE#%i timeseries" % (ii + 1))
        write_split_ts(
            assets.catd[:, :, :, ii, :],
            assets.comptable,
            assets.mmix,
            "e%i" % (ii + 1),
            glsig=glsig,
            assets=assets,
        )


def ctabsel(ctabfile):
    ctlines = open(ctabfile).readlines()
    class_tags = ["#ACC", "#REJ", "#MID", "#IGN"]
    class_dict = {}
    for ii, ll in enumerate(ctlines):
        for kk in class_tags:
            if ll[:4] == kk and ll[4:].strip() != "":
                class_dict[kk] = ll[4:].split("#")[0].split(",")
    return tuple([np.array(class_dict[kk], dtype=int) for kk in class_tags])


###################################################################################################
#                                               Begin Main
###################################################################################################


def me_decompose(assets):
    options = assets.options
    mask = assets.mask
    t2s = assets.t2s
    acc = midk = rej = ignore = []
    if options.mixm is None:
        print("++ Doing ME-PCA and ME-ICA with scikit-learn")
        # import mdp

        if "DEBUG" in sys.argv:
            breakpoint()

        nc, dd = tedpca(options.ste, assets=assets)
        mmix_orig, converge_success = tedica(
            dd, nc, cost=options.initcost, assets=assets
        )

        if not converge_success:
            raise RuntimeError(
                "ICA did not reach any covergence limit. Not producing output."
            )

        np.savetxt("__meica_mix.1D", mmix_orig)
        seldict, comptable, betas, mmix = fitmodels_direct(
            assets.catd,
            mmix_orig,
            mask,
            t2s,
            assets.tes,
            fout=options.fout,
            reindex=True,
            assets=assets,
        )
        np.savetxt("meica_mix.1D", mmix)
        if "GROUP0" in sys.argv:
            group0_flag = True
        else:
            group0_flag = False

        selcomps = SelcompsEncoding(
            seldict=seldict,
            Ne=assets.ne,
            mask=mask,
            t2s=t2s,
            header=assets.head,
            affine=assets.aff,
        )

        acc, rej, midk, ignore = selcomps.fit
        del dd
    else:
        mmix_orig = np.loadtxt("meica_mix.1D")
        # eim = eimask(np.float64(fmask(assets.catd, mask))) == 1
        # eimum = np.array(
        #     np.squeeze(unmask(np.array(eim, dtype=int).prod(1), mask)), dtype=bool
        # )
        seldict, comptable, betas, mmix = fitmodels_direct(
            assets.catd,
            mmix_orig,
            mask,
            t2s,
            assets.tes,
            fout=options.fout,
            assets=assets,
        )
        if options.ctab is None:
            selcomps = SelcompsEncoding(
                seldict=seldict,
                Ne=assets.ne,
                t2s=t2s,
                mask=mask,
                header=assets.head,
                affine=assets.aff,
            )
            acc, rej, midk, ignore = selcomps.fit

    # Do accuracy test if possible
    # selcomps.ground_truth_test()  # type: ignore

    assets.mmix = mmix
    assets.acc = acc
    assets.rej = rej
    assets.midk = midk
    assets.empty = ignore
    assets.comptable = comptable
    assets.seldict = seldict

    if len(acc) == 0:
        print(
            "\n** WARNING! No BOLD components detected!!! Please check data and results!\n"
        )


def tedana_main():
    options, args = tedana_command_line()

    print("-- ME-PCA/ME-ICA Component for ME-ICA %s--" % __version__)

    assets = tedana_process_input(options, args)

    print("++ Computing Mask")
    mask, masksum = makeadmask(assets.catd, min=False, getsum=True)

    print("++ Computing T2* map")
    t2s, s0, t2ss, s0s, t2sG, s0G = t2sadmap(
        assets.catd, mask, assets.tes, masksum, assets=assets
    )

    assets.mask = mask
    assets.t2s = assets.t2sG = t2sG

    # Optimally combine data
    assets.OCcatd = create_memmap("_OCcatd", assets.tsshape)
    assets.OCcatd[:] = optcom(assets.catd, assets.t2s, assets.tes, mask)

    if options.pre_gscontrol:
        assets.gsc_catd = create_memmap("_gsc_catd", assets.catd.shape)
        assets.gsc_catd[:] = gscontrol_raw(OCcatd=assets.OCcatd, assets=assets)

    # TODO: Place daw reducing outer-loop here by counting RuntimeExceptions up to
    # minimum reduced daw ~ 0.5, then do a full fail out
    assets.kdaw_orig = assets.kdaw

    while assets.kdaw > 0.1:
        try:
            print(f"Trying decompositon with kdaw={assets.kdaw}")
            me_decompose(assets)
            break
        except RuntimeError:
            print(f"No ICA solution found at kdaw={assets.kdaw}.")
            assets.kdaw = assets.kdaw * (3.0 / 4.0)
            if assets.kdaw < 0.1:
                raise RuntimeError(
                    "No convergence found across any daws>0.5. Failing. Evaluate data carefully."
                )

    writeresults(assets=assets)

    # breakpoint()

    if options.post_gscontrol:
        glsig = gscontrol_mmix(
            assets=assets,
            OCcatd=assets.OCcatd,
            mmix=assets.mmix,
            mask=assets.mask,
            acc=assets.acc,
        )
    else:
        glsig = None
    assets.glsig = glsig

    writeresults_echoes(glsig=glsig, assets=assets)
    signal_to_dr2s(
        glsig=glsig,
        assets=assets,
    )

    if options.vendor_outputs:
        vendor_outputs(assets=assets)


if __name__ == "__main__":
    tedana_main()
