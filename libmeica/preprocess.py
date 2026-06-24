#!/usr/bin/env python
import subprocess
import sys
# from string import rstrip, split  # python 3 deprecated string import
from optparse import SUPPRESS_HELP, OptionGroup, OptionParser

from packaging import version

from .command_line import _collect_dsinputs

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
# meica.py version %s (c) 2014 Prantik Kundu
# PROCEDURE 1 : Preprocess multi-echo datasets and apply multi-echo ICA based on spatial concatenation
# -Check arguments, input filenames, and filesystem for dependencies
# -Calculation of motion parameters based on images with highest constrast
# -Application of motion correction and coregistration parameters
# -Misc. EPI preprocessing (temporal alignment, smoothing, etc) in appropriate order
# -Compute PCA and ICA in conjuction with TE-dependence analysis
""" % (
    __version__,
    __version__,
)


# Filename parser for NIFTI and AFNI files
def dsprefix(idn):
    # import ipdb; ipdb.set_trace()
    def prefix(datasetname):
        return datasetname.split("+")[0]

    if len(idn.split(".")) != 0:
        if (
            idn.split(".") == "HEAD"
            or idn.split(".")[-1] == "BRIK"
            or idn.split(".")[-2:] == ["BRIK", "gz"]
        ):
            return prefix(idn)
        elif idn.split(".")[-1] == "nii" and not idn.split(".")[-1] == "nii.gz":
            return ".".join(idn.split(".")[:-1])
        elif idn.split(".")[-2:] == ["nii", "gz"]:
            return ".".join(idn.split(".")[:-2])
        else:
            return prefix(idn)
    else:
        return prefix(idn)


def dssuffix(idna):
    suffix = idna.split(dsprefix(idna))[-1]
    # print suffix
    spl_suffix = suffix.split(".")
    # print spl_suffix
    if len(spl_suffix[0]) != 0 and spl_suffix[0][0] == "+":
        return spl_suffix[0]
    else:
        return suffix


def version_checker(cur_ver, ref_ver):
    """
    Checks version in major/minor format of a current version versus a reference version.
    Supports 2 or more level versioning, only 2 levels checked)

    Input:
    cur_ver: float or string of version to check
    ref_ver: float or string of reference version

    Returns:
    bool for pass or fail
    """
    return version.parse(cur_ver) >= version.parse(ref_ver)


# Run dependency check
def dep_check():
    print("++ Checking system for dependencies...")
    fails = 0
    numpy_installed = 0
    scipy_installed = 0
    sklearn_installed = 0
    nibabel_installed = 0
    python_version_ok = 0
    global grayweight_ok
    grayweight_ok = 0
    print(" + Python version: %s" % (".".join([str(v) for v in sys.version_info[0:3]])))
    if sys.version_info < (3, 7):
        print("*+ Python upgrade to Python >= 3.7 & Numpy >= 1.5.x.")
        fails += 1
    else:
        python_version_ok = 1

    try:
        import numpy

        numpy_installed = 1
    except ImportError:
        print(
            "*+ Can't import Numpy! Please check Numpy installation for this Python environment."
        )
        fails += 1

    try:
        import scipy

        scipy_installed = 1
    except ImportError:
        print(
            "*+ Can't import Scipy! Please check Scipy installation for this Python environment."
        )
        fails += 1

    try:
        import sklearn

        sklearn_installed = 1
    except ImportError:
        print(
            "*+ Can't import scikit-learn! Please install scikit-learn version >0.15.0 for this Python environment."
        )
        fails += 1

    try:
        import nibabel

        nibabel_installed = 1
    except ImportError:
        print("*+ Can't import Nibabel! Please install for this Python environment.")
        fails += 1

    if numpy_installed:
        print(" + Numpy version: %s" % (numpy.__version__))  # type: ignore
        if version_checker(numpy.__version__, "1.5") is False:  # type: ignore
            fails += 1
            print("*+ Numpy version is too old! Please upgrade to Numpy >=1.5.x!")
        # Check BLAS installation (python 3)
        # blas_check = True in ["blas" in ff for ff in dir(numpy.__config__)]
        # if not blas_check:
        #     fails += 1
        #     print("*+ Numpy is not linked to BLAS! Please check Numpy installation.")
    if scipy_installed:
        print(" + Scipy version: %s" % (scipy.__version__))  # type: ignore
        if version_checker(scipy.__version__, "0.11") is False:  # type: ignore
            fails += 1
            print("*+ Scipy version is too old! Please upgrade to Scipy >=0.11.x!")
    if sklearn_installed:
        print(" + scikit-learn version: %s" % (sklearn.__version__))  # type: ignore
        if version_checker(sklearn.__version__, "0.15") is False:  # type: ignore
            fails += 1
            print(
                "*+ scikit-learn version is too old! Please upgrade to Scipy >=0.15.x!"
            )
    if nibabel_installed:
        print(" + nibabel version: %s" % (nibabel.__version__))  # type: ignore
        if version_checker(nibabel.__version__, "2.4.0") is False:  # type: ignore
            fails += 1
            print("*+ nibabel version is too old! Please 'pip install nibabel==2.4.0")
    afnicheck = subprocess.getstatusoutput("3dinfo")
    afnisegcheck = subprocess.getstatusoutput("3dSeg -help")
    if afnicheck[0] != 0:
        print("*+ Can't run AFNI binaries. Make sure AFNI is on the path!")
        fails += 1
    elif (
        not afnicheck[1].__contains__("Alternate Alternative Usage")
        and len(afnisegcheck[1]) < 1000
    ):
        print(
            "*+ This seems like an old version of AFNI. Please upgrade to latest version of AFNI."
        )
        fails += 1
    if afnisegcheck[0] == 0 and len(afnisegcheck[1]) >= 1000:
        print(
            " + Using AFNI 3dSeg for gray matter weighted anatomical-functional coregistration"
        )
        grayweight_ok = 1
    if grayweight_ok == 0:
        print(
            "*+ WARNING: AFNI 3dSeg not available for gray matter weighted coregistration. See README."
        )
    if fails == 0:
        print(" + Dependencies OK.")
    else:
        print("*+ EXITING. Please see error messages.")
        sys.exit(1)


def argsort(seq):
    return sorted(range(len(seq)), key=seq.__getitem__)


def logcomment(comment, level=3, buffer=None, global_sl=None):
    if global_sl is not None:
        global sl
        sl = global_sl
    if buffer is None:
        buffer = sl
    out_cmd = "echo"
    majmark = "\n"
    leading = "--------"
    if level == 5:
        # A regular comment with an echo
        majmark = ""
        leading = ""
    if level == 4:
        # A regular comment
        majmark = "#"
        leading = ""
        out_cmd = ""
    if level == 3:
        # A simple echo
        majmark = ""
    if level == 2:
        # A section delimiter
        leading = "\n--------\n"
    if level == 1:
        # An echo with AFNI-style emphasis
        leading = "+* "
        buffer.append("""echo "\n++++++++++++++++++++++++" """)
        majmark = ""
    buffer.append(f'{majmark}{out_cmd} "{leading} {comment}"')


# Configure options and help dialog
parser = OptionParser()
parser.add_option(
    "-e",
    "",
    dest="tes",
    action="callback",
    callback=_collect_dsinputs,
    type="string",
    default=["bids"],
    help="ex: -e bids | -e 14.5 38.5 62.5  Echo times from BIDS header or in ms",
)
parser.add_option(
    "-d",
    "",
    dest="dsinputs",
    action="callback",
    callback=_collect_dsinputs,
    type="string",
    help="ex: -d RESTe1.nii.gz RESTe2.nii.gz RESTe3.nii.gz",
    default=[],
)
parser.add_option(
    "-a",
    "",
    dest="anat",
    help="ex: -a mprage.nii.gz  Anatomical dataset (optional)",
    default="",
)
parser.add_option(
    "-b",
    "",
    dest="basetime",
    help="ex: -b 10s OR -b 10v  Time to steady-state equilibration in seconds(s) or volumes(v). Default 0. ",
    default="0",
)
parser.add_option(
    "",
    "--MNI",
    dest="mni",
    action="store_true",
    help="Warp to MNI standard space using high-resolution template",
    default=False,
)
extopts = OptionGroup(parser, "Additional processing options")
extopts.add_option(
    "",
    "--qwarp",
    dest="qwarp",
    action="store_true",
    help="Nonlinear warp to standard space using QWarp (use --MNI or --space)",
    default=False,
)
extopts.add_option(
    "",
    "--native",
    dest="native",
    action="store_true",
    help="Output native space results in addition to standard space results",
    default=False,
)
extopts.add_option(
    "",
    "--space",
    dest="space",
    help="Path to specific standard space template for affine anatomical normalization",
    default=False,
)
extopts.add_option(
    "",
    "--fres",
    dest="fres",
    help="Specify functional voxel dim. in mm (iso.) for resampling during preprocessing. Default none. ex: --fres=2.5",
    default=False,
)
extopts.add_option(
    "",
    "--ted_spms",
    dest="ted_spms",
    help="Export TE-dependence SPMs (warning: many files)",
    action="store_true",
    default=False,
)
extopts.add_option(
    "",
    "--ted_seed",
    dest="ted_seed",
    help="Seed for determinisitc / pseudo-random FastICA. Default fully random. ex: --ted_seed=123",
    default=False,
)
extopts.add_option(
    "",
    "--no_skullstrip",
    action="store_true",
    dest="no_skullstrip",
    help="Anatomical is already intensity-normalized and skull-stripped (if -a provided)",
    default=False,
)
extopts.add_option(
    "",
    "--no_despike",
    action="store_true",
    dest="no_despike",
    help="Do not de-spike functional data. Default is to de-spike, recommended.",
    default=False,
)
extopts.add_option(
    "",
    "--axialize",
    action="store_true",
    dest="axialize",
    help="Do not re-write dataset in axial-first order. Default is to axialize, recommended.",
    default=False,
)
extopts.add_option(
    "",
    "--mask_mode",
    dest="mask_mode",
    help="Mask functional with help from anatomical or standard space images: use 'anat' or 'template' ",
    default="func",
)
extopts.add_option(
    "",
    "--coreg_mode",
    dest="coreg_mode",
    help="Coregistration with Local Pearson and T2* weights (default), or use align_epi_anat.py (edge method): use 'lp-t2s' or 'aea'",
    default="lp-t2s",
)
# parser.add_option('',"--strict",dest='strict',action='store_true',help="Use strict component selection, suitable to large-voxel high-SNR data",default=False)
parser.add_option(
    "",
    "--strict",
    dest="strict",
    action="store_true",
    help=SUPPRESS_HELP,
    default=False,
)
extopts.add_option(
    "--smooth",
    dest="FWHM",
    help="Data FWHM smoothing (3dBlurInMask). Default off. ex: --smooth 3mm ",
    default="0mm",
)
extopts.add_option(
    "",
    "--align_base",
    dest="align_base",
    help="Explicitly specify base dataset for volume registration",
    default="",
)
extopts.add_option(
    "",
    "--TR",
    dest="TR",
    help="The TR. Default read from input dataset header",
    default="",
)
extopts.add_option(
    "",
    "--tpattern",
    dest="tpattern",
    help="Slice timing (i.e. alt+z, see 3dTshift -help). Default from header. (N.B. This is important!)",
    default="",
)
extopts.add_option(
    "",
    "--align_args",
    dest="align_args",
    help="Additional arguments to anatomical-functional co-registration routine",
    default="",
)
extopts.add_option(
    "",
    "--ted_args",
    dest="ted_args",
    help="Additional arguments to TE-dependence analysis routine",
    default="",
)
extopts.add_option(
    "",
    "--select_only",
    dest="select_only",
    action="store_true",
    help=SUPPRESS_HELP,
    default=False,
)
extopts.add_option(
    "",
    "--tedica_only",
    dest="tedica_only",
    action="store_true",
    help=SUPPRESS_HELP,
    default=False,
)
extopts.add_option(
    "",
    "--export_only",
    dest="export_only",
    action="store_true",
    help=SUPPRESS_HELP,
    default=False,
)
extopts.add_option("", "--daw", dest="daw", help=SUPPRESS_HELP, default="100")
extopts.add_option(
    "", "--tlrc", dest="space", help=SUPPRESS_HELP, default=False
)  # For backwards compat. with existing scripts
extopts.add_option("", "--highpass", dest="highpass", help=SUPPRESS_HELP, default=0.0)
extopts.add_option("", "--detrend", dest="detrend", help=SUPPRESS_HELP, default=0.0)
extopts.add_option(
    "", "--initcost", dest="initcost", help=SUPPRESS_HELP, default="tanh"
)
extopts.add_option(
    "", "--finalcost", dest="finalcost", help=SUPPRESS_HELP, default="tanh"
)
extopts.add_option(
    "", "--sourceTEs", dest="sourceTEs", help=SUPPRESS_HELP, default="-1"
)
parser.add_option_group(extopts)
runopts = OptionGroup(parser, "Run options")
runopts.add_option(
    "",
    "--work_root",
    dest="work_root",
    help="Directory to write meicadir during processing",
    default="",
)
runopts.add_option(
    "",
    "--no_trim_work",
    dest="no_trim_work",
    help="Retain all intermediate files made during preprocessing.",
    action="store_true",
    default=False,
)
runopts.add_option(
    "",
    "--prefix",
    dest="prefix",
    help="Prefix for final ME-ICA output datasets.",
    default="",
)
runopts.add_option(
    "-j",
    "--cpus",
    dest="cpus",
    help="Maximum number of CPUs (OpenMP threads) to use. Default 4.",
    default="4",
)
runopts.add_option(
    "",
    "--no_para_echo",
    dest="no_para_echo",
    action="store_true",
    help="Parallel processing of echo datasets.",
    default=False,
)
runopts.add_option(
    "", "--label", dest="label", help="Label to tag ME-ICA analysis folder.", default=""
)
runopts.add_option(
    "",
    "--test_proc",
    action="store_true",
    dest="test_proc",
    help="Align and preprocess 1 dataset then exit, for testing",
    default=False,
)
runopts.add_option(
    "",
    "--script_only",
    action="store_true",
    dest="script_only",
    help="Generate script only, then exit",
    default=0,
)
runopts.add_option(
    "",
    "--pp_only",
    action="store_true",
    dest="pp_only",
    help="Preprocess only, then exit.",
    default=False,
)
runopts.add_option(
    "",
    "--keep_int",
    action="store_true",
    dest="keep_int",
    help="Keep preprocessing intermediates. Default delete.",
    default=False,
)
runopts.add_option(
    "",
    "--skip_check",
    action="store_true",
    dest="skip_check",
    help="Skip dependency checks during initialization.",
    default=False,
)
runopts.add_option(
    "",
    "--RESUME",
    dest="resume",
    action="store_true",
    help="Attempt to resume from normalization step onwards, overwriting existing",
)
runopts.add_option(
    "",
    "--OVERWRITE",
    dest="overwrite",
    action="store_true",
    help="If meica.xyz directory exists, overwrite. ",
    default=False,
)
runopts.add_option(
    "",
    "--OVERWRITE_NOWAIT",
    dest="overwrite_nowait",
    action="store_true",
    help="",
    default=False,
)
parser.add_option_group(runopts)
(options, args) = parser.parse_args()
