import atexit
import os
import sys
from optparse import OptionParser
from pathlib import Path
from subprocess import CalledProcessError
from types import SimpleNamespace

import nibabel as nib
import numpy as np

from .utils.volume import cat2echos


def _unlink(p):
    p = Path(p)
    p.unlink()


def tedana_command_line():
    parser = OptionParser()
    parser.add_option(
        "-d",
        "--orig_data",
        dest="data",
        help="Spatially Concatenated Multi-Echo Dataset",
        default=None,
    )
    parser.add_option(
        "-e", "--TEs", dest="tes", help="Echo times (in ms) ex: 15,39,63", default=None
    )
    parser.add_option(
        "",
        "--mix",
        dest="mixm",
        help="Mixing matrix. If not provided, ME-PCA & ME-ICA (MDP) is done.",
        default=None,
    )
    parser.add_option(
        "",
        "--ctab",
        dest="ctab",
        help="Component table extract pre-computed classifications from.",
        default=None,
    )
    parser.add_option(
        "",
        "--T2_gray",
        dest="T2_gray",
        help="T2 of gray matter, usually 80 at 3T and 50 at 7T",
        default=80,
        type="float",
    )
    parser.add_option(
        "",
        "--manacc",
        dest="manacc",
        help="Comma separated list of manually accepted components",
        default=None,
    )
    parser.add_option(
        "",
        "--strict",
        dest="strict",
        action="store_true",
        help="Ignore low-variance ambiguous components",
        default=False,
    )
    # parser.add_option('',"--wav",dest='wav',help="Perform wavelet PCA, default False",default=False)
    parser.add_option(
        "",
        "--pre_gscontrol",
        dest="pre_gscontrol",
        action="store_true",
        help="Control spatial global signal in inputs (i.e., for 7T)",
        default=False,
    )
    parser.add_option(
        "",
        "--post_gscontrol",
        dest="post_gscontrol",
        action="store_true",
        help="Control spatial global signal in ouputs (i.e., for 3T)",
        default=False,
    )
    parser.add_option(
        "",
        "--kdaw",
        dest="kdaw",
        help="Dimensionality augmentation weight (Kappa). Default 10. -1 for low-dimensional ICA",
        default=10.0,
    )
    parser.add_option(
        "",
        "--rdaw",
        dest="rdaw",
        help="Dimensionality augmentation weight (Rho). Default 1. -1 for low-dimensional ICA",
        default=1.0,
    )
    parser.add_option(
        "",
        "--conv",
        dest="conv",
        help="Convergence limit. Default 2.5e-5",
        default="2.5e-5",
    )
    parser.add_option(
        "",
        "--sourceTEs",
        dest="ste",
        help="Source TEs for models. ex: -ste 2,3 ; -ste 0 for all, -1 for opt. com. Default -1.",
        default=0 - 1,
    )
    parser.add_option(
        "",
        "--combmode",
        dest="combmode",
        help="Combination scheme for TEs: t2s (Posse 1999, default),ste(Poser)",
        default="t2s",
    )
    parser.add_option(
        "",
        "--denoiseTEs",
        dest="dne",
        action="store_true",
        help="Denoise each TE dataset separately",
        default=False,
    )
    parser.add_option(
        "",
        "--initcost",
        dest="initcost",
        help="Initial cost func. for ICA: pow3,tanh(default),gaus,skew",
        default="tanh",
    )
    parser.add_option(
        "",
        "--finalcost",
        dest="finalcost",
        help="Final cost func, same opts. as initial",
        default="tanh",
    )
    parser.add_option(
        "",
        "--stabilize",
        dest="stabilize",
        action="store_true",
        help="Stabilize convergence by reducing dimensionality, for low quality data",
        default=False,
    )
    parser.add_option(
        "",
        "--fout",
        dest="fout",
        help="Output accepted TE-dependence SPMs",
        action="store_true",
        default=False,
    )
    parser.add_option(
        "",
        "--fout_all",  # not implemented!!
        dest="fout_all",
        help="Output all TE-dependence Kappa/Rho SPMs",
        action="store_true",
        default=False,
    )
    parser.add_option(
        "",
        "--filecsdata",
        dest="filecsdata",
        help="Save component selection data",
        action="store_true",
        default=False,
    )
    parser.add_option(
        "", "--label", dest="label", help="Label for output directory.", default=None
    )

    (options, args) = parser.parse_args()
    args = list(set(args) - set(["DEBUG"]))

    return options, args


def tedana_process_input(options, args):
    if options.tes is None or options.data is None:
        print("*+ Need at least data and TEs, use -h for help.")
        sys.exit()

    print("++ Loading Data")
    tes = np.fromstring(options.tes, sep=",", dtype=np.float32)
    ne = tes.shape[0]
    catim = nib.load(options.data)  # type: ignore
    head = catim.header
    head.extensions = []  # type: ignore
    head.set_sform(head.get_sform(), code=1)  # type: ignore
    aff = catim.affine  # type: ignore
    memmap_path = "_catd.memmap"
    catim_shape = catim.shape  # type: ignore
    catd_shape = catim_shape[:2] + (catim_shape[2] // ne,) + (ne,) + (catim_shape[3],)
    catd = np.memmap(memmap_path, dtype=np.float32, mode="w+", shape=catd_shape)  # type: ignore
    atexit.register(_unlink, memmap_path)
    # import pudb; pudb.set_trace()
    catd[:] = cat2echos(catim.get_fdata(), ne)  # changed here too  # type: ignore
    catd.flush()
    nx, ny, nz, Ne, nt = catd.shape
    mu = catd.mean(axis=-1)
    sig = catd.std(axis=-1)

    """Parse options, prepare output directory"""
    # if options.fout:
    #     options.fout = aff
    # else:
    #     options.fout = None
    kdaw = float(options.kdaw)
    rdaw = float(options.rdaw)
    if options.label is not None:
        dirname = "%s" % ".".join(["TED", options.label])
    else:
        dirname = "TED"
    os.system("mkdir %s" % dirname)
    if options.mixm is not None:
        try:
            os.system(
                "cp %s %s/meica_mix.1D; cp %s %s/%s"
                % (
                    options.mixm,
                    dirname,
                    options.mixm,
                    dirname,
                    os.path.basename(options.mixm),
                )
            )
        except CalledProcessError:
            pass
    if options.ctab is not None:
        try:
            os.system(
                "cp %s %s/comp_table.txt; cp %s %s/%s"
                % (
                    options.mixm,
                    dirname,
                    options.mixm,
                    dirname,
                    os.path.basename(options.mixm),
                )
            )
        except CalledProcessError:
            pass

    os.chdir(dirname)
    assets = SimpleNamespace()
    assets.nt = nt
    assets.ne = ne
    assets.Ne = ne
    assets.head = head
    assets.aff = aff
    assets.catd = catd
    assets.mu = mu
    assets.sig = sig
    assets.dirname = dirname
    assets.tes = tes
    assets.kdaw = kdaw
    assets.rdaw = rdaw
    assets.args = args
    assets.tsshape = (nx, ny, nz, nt)
    assets.options = options
    return assets
