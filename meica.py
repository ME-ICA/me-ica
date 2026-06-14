#!/usr/bin/env python
import json
import os.path
import random
import re
import subprocess
import sys

# from string import rstrip, split  # python 3 deprecated string import
from optparse import SUPPRESS_HELP, OptionGroup, OptionParser
from os import chdir, getcwd, mkdir, popen, system
from re import split as resplit

from libmeica.preprocess import (
    argsort,
    dep_check,
    dsprefix,
    dssuffix,
    logcomment,
    options,
)

VENDOR_OPT = None
VENDOR_ARGS = None
try:
    from libmeica.vendor.output import vendor_exports, vendortedargs

    VENDOR_OPT = "--vendor_output"
    VENDOR_ARGS = vendortedargs()
except ImportError:
    pass

global sl
sl = []  # Script command list

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


# Welcome line
print(
    """\n-- Multi-Echo Independent Components Analysis (ME-ICA) %s --

Please cite:
Kundu, P., Brenowitz, N.D., Voon, V., Worbe, Y., Vertes, P.E., Inati, S.J., Saad, Z.S.,
Bandettini, P.A. & Bullmore, E.T. Integrated strategy for improving functional
connectivity mapping using multiecho fMRI. PNAS (2013).

Kundu, P., Inati, S.J., Evans, J.W., Luh, W.M. & Bandettini, P.A. Differentiating
BOLD and non-BOLD signals in fMRI time series using multi-echo EPI. NeuroImage (2011).
"""
    % (__version__)
)


def getdsname(e_ii, *, prefixonly=False, wsrc=False):
    global STARTDIR
    srcdir = STARTDIR
    if wsrc:
        _srcdir = SRC_FOLDERS[e_ii]
        if len(_srcdir) > 0:
            if _srcdir[0] == "/":
                srcdir = _srcdir
            else:
                srcdir = os.path.join(STARTDIR, _srcdir)
    if SHORTHAND_DSIN:
        dsname = "%s%s%s%s" % (PREFIX, DATASETS[e_ii], TRAILING, ISF)
    else:
        dsname = DATASETS_IN[e_ii]
    if prefixonly:
        if wsrc:
            return os.path.join(srcdir, dsprefix(dsname))
        else:
            return dsprefix(dsname)
    else:
        if wsrc:
            return os.path.join(srcdir, dsname)
        else:
            return dsname


# Parse dataset input names
if options.dsinputs == "" or options.TR == 0:
    if not options.skip_check:
        dep_check()
    print("*+ Need at least dataset inputs and TE. Try meica.py -h")
    sys.exit(1)
if os.path.abspath(os.path.curdir).__contains__("meica."):
    print(
        "*+ You are inside a ME-ICA directory! Please leave this directory and rerun."
    )
    sys.exit(1)


def parse_bids_tes(datasets):
    tes = []
    for _di, _dn in enumerate(datasets):
        _dsin = getdsname(_di, wsrc=True)
        try:
            json_file = f"{getdsname(_di, prefixonly=True, wsrc=True)}.json"
            with open(json_file) as bids_ifh:
                _bids_in = json.load(bids_ifh)
                _echo_time = float(_bids_in["EchoTime"])
                if _echo_time < 1:
                    _echo_time = _echo_time * 1000
                tes.append(str(_echo_time))
        except Exception as e:
            print(e)
            print("*+ Could not load TEs from BIDS!")
            sys.exit(10)
    return tes


if len(sys.argv) == 1:
    dep_check()
    exit()

outprefix = options.prefix

if options.anat == "" and (options.mni or options.space):
    print("Need to specify anatomical for standard space normalization.")
    sys.exit(1)


if isinstance(options.dsinputs, list) and len(options.dsinputs) > 1:
    DATASETS_IN = options.dsinputs
    SRC_FOLDERS = [os.path.dirname(_d) for _d in DATASETS_IN]
    DATASETS_IN = [os.path.basename(_d) for _d in DATASETS_IN]
    DATASETS = [str(vv + 1) for vv in range(len(options.dsinputs))]
    PREFIX = dsprefix(DATASETS_IN[0])
    ISF = dssuffix(DATASETS_IN[0])
    if ".nii" in ISF:
        ISF = ".nii"
    TRAILING = ""
    SETNAME = PREFIX + options.label
    options.dsinputs = ",".join(DATASETS_IN)
    SHORTHAND_DSIN = False
else:
    if isinstance(options.dsinputs, list):
        options.dsinputs = options.dsinputs[0]
    if "[" in options.dsinputs:
        SHORTHAND_DSIN = True
        dsinputs = dsprefix(options.dsinputs)
        PREFIX = resplit(r"[\[\],]", dsinputs)[0]
        DATASETS = resplit(r"[\[\],]", dsinputs)[1:-1]
        TRAILING = resplit(r"[\]+]", dsinputs)[-1]
        ISF = dssuffix(options.dsinputs)
        SETNAME = PREFIX + "".join(DATASETS) + TRAILING + options.label
        SRC_FOLDERS = [os.path.dirname(_d) for _d in DATASETS]
        DATASETS_IN = DATASETS = [os.path.basename(_d) for _d in DATASETS]
    else:
        # Parse longhand input file specificiation
        SHORTHAND_DSIN = False
        DATASETS_IN = options.dsinputs.split(",")
        assert isinstance(DATASETS_IN, list)
        SRC_FOLDERS = [os.path.dirname(_d) for _d in DATASETS_IN]
        DATASETS_IN = [os.path.basename(_d) for _d in DATASETS_IN]
        DATASETS = [str(vv + 1) for vv in range(len(DATASETS_IN))]
        PREFIX = dsprefix(DATASETS_IN[0])
        ISF = dssuffix(DATASETS_IN[0])
        if ".nii" in ISF:
            ISF = ".nii"
        TRAILING = ""
        SETNAME = PREFIX + options.label

# Prepare script
STARTDIR = str.rstrip(popen("pwd").readlines()[0])

# Parse shorthand input file specification and TEs
if isinstance(options.tes, list) and len(options.tes) > 1:
    tes = options.tes[:]
else:
    if isinstance(options.tes, list):
        options.tes = options.tes[0]
    if "," in options.tes:
        tes = str(options.tes).split(",")
    else:
        tes = parse_bids_tes(DATASETS)
tes_float = [float(_te) for _te in tes]
tes_order = argsort(tes_float)
DATASETS_IN = [DATASETS_IN[_tei] for _tei in tes_order]
tes = [tes_float[_tei] for _tei in tes_order]
options.tes = ",".join(["%0.2f" % _te for _te in tes])

# import pudb
#
# pudb.set_trace()


if not SHORTHAND_DSIN and len(DATASETS) != len(DATASETS_IN):  # type: ignore
    print(
        "*+ Can't understand dataset specification. Try double quotes around -d argument."
    )
    sys.exit(1)

if len(options.tes.split(",")) != len(DATASETS):
    print(
        "*+ Number of TEs and input datasets must be equal and matched in order. Or try double quotes around -d argument."
    )
    sys.exit(1)


# Determine echo processing parallelism
para_mode = not options.no_para_echo
para_echo_indicator = ""
omp_threads_line = f"OMP_NUM_THREADS={options.cpus}"
omp_threads_line_dspk = omp_threads_line
if para_mode:
    para_echo_indicator = "&"
    omp_threads_line = (
        f"OMP_NUM_THREADS={max(1, round(int(options.cpus) / len(options.tes)))}"
    )
    omp_threads_line_dspk = (
        f"OMP_NUM_THREADS={min(4, max(1, round(int(options.cpus) / len(options.tes))))}"
    )
omp_threads_line_dspk_max = f"OMP_NUM_THREADS={min(int(options.cpus), 4)}"
omp_threads_half = max(round(float(options.cpus) / 2), 1)


# Prepare alternate work dir location during processing
linked_workdir = False
if options.work_root != "":
    work_root = os.path.abspath(options.work_root)
    linked_workdir = True
    workdir_hash = random.getrandbits(128)
    workdir_full_path = f"{work_root}/meica.{SETNAME}.{workdir_hash}"
else:
    work_root = os.path.abspath(".")
    workdir_full_path = f"{work_root}/meica.%s" % (SETNAME)

finaldir_full_path = f"{STARTDIR}/meica.%s" % (SETNAME)

meica_src_dir = os.path.dirname(os.path.abspath(os.path.expanduser(sys.argv[0])))
headsl = []  # Header lines and command list
runcmd = (
    " ".join(sys.argv)
    .replace(options.dsinputs, r'"%s"' % options.dsinputs)
    .replace('"', r'"')
)
headsl.append("#" + runcmd)
headsl.append(welcome_block)
osf = ".nii.gz"  # Using NIFTI outputs

logcomment("Parameters:", buffer=headsl, level=2)
logcomment(f"Datasets: {options.dsinputs}", buffer=headsl, level=5)
logcomment(f"TEs: {options.tes}", buffer=headsl, level=5)
# Check if input files exist

notfound = 0
for ds_ii in range(len(DATASETS)):
    if subprocess.getstatusoutput("3dinfo %s" % (getdsname(ds_ii, wsrc=True)))[0] != 0:
        print("*+ Can't find/load dataset %s !" % (getdsname(ds_ii, wsrc=True)))
        notfound += 1
if (
    options.anat != ""
    and subprocess.getstatusoutput("3dinfo %s" % (options.anat))[0] != 0
):
    print("*+ Can't find/load anatomical dataset %s !" % (options.anat))
    notfound += 1
if notfound != 0:
    print("++ EXITING. Check dataset names.")
    sys.exit(1)

# Check dependencies
grayweight_ok = 0
if not options.skip_check:
    dep_check()
    print("++ Continuing with preprocessing.")
else:
    print("*+ Skipping dependency checks.")
    grayweight_ok = 1

# Parse timing arguments
if options.TR != "":
    tr = float(options.TR)
else:
    tr = float(
        os.popen("3dinfo -tr %s" % (getdsname(0, wsrc=True))).readlines()[0].strip()
    )
    options.TR = str(tr)
if "v" in str(options.basetime):
    basebrik = int(options.basetime.strip("v"))
else:
    timetoclip = 0
    timetoclip = float(options.basetime.strip("s"))
    basebrik = int(round(timetoclip / tr))

# Misc. command parsing
if options.mni:
    options.space = "MNI_caez_N27+tlrc"
if options.qwarp and (options.anat == "" or not options.space):
    print(
        "*+ Can't specify Qwarp nonlinear coregistration without anatomical and SPACE template!"
    )
    sys.exit(1)

if not options.mask_mode in ["func", "anat", "template"]:
    print("*+ Mask mode option '%s' is not recognized!" % options.mask_mode)
    sys.exit(1)
if options.mask_mode == "" and options.space:
    options.mask_mode = "template"

# Parse alignment options
if options.coreg_mode == "aea":
    options.t2salign = False
elif "lp" in options.coreg_mode:
    options.t2salign = True
align_base = basebrik
align_interp = "cubic"
align_interp_final = "wsinc5"
oblique_epi_read = 0
oblique_anat_read = 0
zeropad_opts = " -I %s -S %s -A %s -P %s -L %s -R %s " % (tuple([1] * 6))
if options.anat != "":
    oblique_anat_read = int(
        os.popen("3dinfo -is_oblique %s" % (options.anat)).readlines()[0].strip()
    )
    epicm = [
        float(coord)
        for coord in os.popen("3dCM %s" % (getdsname(0, wsrc=True)))
        .readlines()[0]
        .strip()
        .split()
    ]
    anatcm = [
        float(coord)
        for coord in os.popen("3dCM %s" % (options.anat)).readlines()[0].strip().split()
    ]
    maxvoxsz = float(
        os.popen("3dinfo -dk %s" % (getdsname(0, wsrc=True))).readlines()[0].strip()
    )
    deltas = [
        abs(epicm[0] - anatcm[0]),
        abs(epicm[1] - anatcm[1]),
        abs(epicm[2] - anatcm[2]),
    ]
    cmdist = 20 + sum([dd**2.0 for dd in deltas]) ** 0.5
    cmdif = max(
        abs(epicm[0] - anatcm[0]), abs(epicm[1] - anatcm[1]), abs(epicm[2] - anatcm[2])
    )
    addslabs = abs(int(cmdif / maxvoxsz)) + 10
    zeropad_opts = " -I %s -S %s -A %s -P %s -L %s -R %s " % (tuple([addslabs] * 6))
oblique_epi_read = int(
    os.popen("3dinfo -is_oblique %s" % (getdsname(0, wsrc=True))).readlines()[0].strip()
)
if oblique_epi_read or oblique_anat_read:
    oblique_mode = True
    headsl.append("echo Oblique data detected.")
else:
    oblique_mode = False
if options.fres:
    if options.qwarp:
        qwfres = "-dxyz %s" % options.fres
    alfres = "-mast_dxyz %s" % options.fres
else:
    if options.qwarp:
        qwfres = "-dxyz ${voxsize}"  # See section called "Preparing functional masking for this ME-EPI run"
    alfres = "-mast_dxyz ${voxsize}"
if options.anat == "" and options.mask_mode != "func":
    print("*+ Can't do anatomical-based functional masking without an anatomical!")
    sys.exit(1)
if options.anat and options.space and options.qwarp:
    valid_qwarp_mode = True
else:
    valid_qwarp_mode = False

# Detect if current AFNI has old 3dNwarpApply
if (
    " -affter aaa  = *** THIS OPTION IS NO LONGER AVAILABLE"
    in subprocess.getstatusoutput("3dNwarpApply -help")[1]
):
    old_qwarp = False
else:
    old_qwarp = True

# Detect AFNI direcotry
afnidir = os.path.dirname(os.popen("which 3dSkullStrip").readlines()[0])

# Prepare script and enter MEICA directory
logcomment("Set up script run environment", level=1, global_sl=sl)
headsl.append("set -e")
headsl.append("export OMP_NUM_THREADS=%s" % (options.cpus))
headsl.append("export MKL_NUM_THREADS=%s" % (options.cpus))
headsl.append("export DYLD_FALLBACK_LIBRARY_PATH=%s" % (afnidir))
headsl.append("export AFNI_3dDespike_NEW=YES")

if options.export_only or options.select_only or options.tedica_only:
    headsl.append(f"voxdims=$(cat meica.{SETNAME}/voxdims.1D)")
    headsl.append(f"voxsize=$(cat meica.{SETNAME}/voxsize.1D)")

initial_run = (
    not options.resume
    and not options.tedica_only
    and not options.select_only
    and not options.export_only
)

# Clear workdir
if linked_workdir:
    headsl.append(f"rm -rf {workdir_full_path}")

if initial_run:
    if options.overwrite or options.overwrite_nowait:  # if we're overwriting
        if options.overwrite and options.overwrite_nowait:  # if both overwrites,
            print("Specify either --OVERWRITE or --OVERWRITE_NOWAIT; not both.")
            sys.exit(1)  # bail
        else:
            if (
                options.overwrite and not options.overwrite_nowait
            ):  # if overwrite with wait
                logcomment(
                    "WAITING FOR 5 SECONDS BEFORE OVERWRITING. Ctrl-C TO ABORT.",
                    1,
                    buffer=headsl,
                )
                headsl.append("sleep 5")
            headsl.append(f"rm -rf {finaldir_full_path}")
    else:
        headsl.append(
            f"if [[ -e {finaldir_full_path} || -L {finaldir_full_path} ]]; then echo ME-ICA directory exists, exiting; exit; fi"
        )
    headsl.append(f"mkdir -p {workdir_full_path}")
    if linked_workdir:
        headsl.append(f"ln -s {workdir_full_path} {finaldir_full_path}")

if options.resume:
    headsl.append(
        "if [ ! -e meica.%s/_meica.orig.sh ]; then mv `ls meica.%s/_meica*sh` meica.%s/_meica.orig.sh; fi"
        % (SETNAME, SETNAME, SETNAME)
    )
if not options.tedica_only and not options.select_only:
    headsl.append("cp _meica_%s.sh meica.%s/" % (SETNAME, SETNAME))

headsl.append("cd meica.%s" % SETNAME)
thecwd = "%s/meica.%s" % (getcwd(), SETNAME)

ica_datasets = sorted(DATASETS)

# Parse anatomical processing options, process anatomical
ANAT_DIR = None
if options.anat != "":
    sl.append(f"anat_proc_st1()" + "{")
    logcomment(
        "Deoblique, unifize, skullstrip, and/or autobox anatomical, in starting directory (may take a little while)",
        level=1,
    )
    nsmprage = options.anat
    ANAT_DIR = os.path.dirname(nsmprage) if "/" in nsmprage else STARTDIR
    if ANAT_DIR[0] != "/" and ANAT_DIR != STARTDIR:
        ANAT_DIR = os.path.join(STARTDIR, ANAT_DIR)
    nsmprage = os.path.basename(options.anat)
    anatprefix = dsprefix(nsmprage)
    pathanatprefix = "%s/%s" % (ANAT_DIR, anatprefix)
    if oblique_mode:
        sl.append(
            "if [ ! -e %s_do.nii.gz ]; then 3dWarp -wsinc5 -overwrite -prefix %s_do.nii.gz -deoblique %s/%s; fi"
            % (pathanatprefix, pathanatprefix, ANAT_DIR, nsmprage)
        )
        nsmprage = "%s_do.nii.gz" % (anatprefix)
    if not options.no_skullstrip:
        sl.append(
            "if [ ! -e %s_ns.nii.gz ]; then 3dUnifize -overwrite -prefix %s_u.nii.gz %s/%s; 3dSkullStrip  -shrink_fac_bot_lim 0.3 -orig_vol -overwrite -prefix %s_wbns.nii.gz -input %s_u.nii.gz; 3dAutobox -overwrite -prefix %s_ns.nii.gz %s_wbns.nii.gz; fi"
            % (
                pathanatprefix,
                pathanatprefix,
                ANAT_DIR,
                nsmprage,
                pathanatprefix,
                pathanatprefix,
                pathanatprefix,
                pathanatprefix,
            )
        )
        nsmprage = "%s_ns.nii.gz" % (anatprefix)
    sl.append("}")

    sl.append(f"anat_proc_st1 {para_echo_indicator}")

# Copy in functional datasets as NIFTI (if not in NIFTI already), calculate rigid body alignment
# import ipdb; ipdb.set_trace()
vrbase = getdsname(0, prefixonly=True)
logcomment("Copy in functional datasets, reset NIFTI tags as needed", level=1)
for e_ii in range(len(DATASETS)):
    sl.append(
        "3dcalc -a %s -expr 'a' -prefix %s.nii"
        % (getdsname(e_ii, wsrc=True), getdsname(e_ii, prefixonly=True))
    )
    if ".nii" in ISF:
        sl.append(
            "nifti_tool -mod_hdr -mod_field sform_code 1 -mod_field qform_code 1 -infiles ./%s.nii -overwrite"
            % (getdsname(e_ii, prefixonly=True))
        )
ISF = ".nii"

logcomment(
    "Calculate and save motion and obliquity parameters, despiking first if not disabled, and separately save and mask the base volume",
    level=1,
)
# Determine input to volume registration
vrinput = "./%s%s" % (vrbase, ISF)
vrAinput = "./%s_vrA%s" % (vrbase, osf)
# Compute obliquity matrix
if oblique_mode:
    # if options.anat != "":
    if options.anat == "":
        sl.append(
            "3dWarp -wsinc5 -overwrite -prefix %s -deoblique %s" % (vrAinput, vrinput)
        )
# Despike and axialize
if not options.no_despike:
    if options.anat != "":
        _despike_input = vrinput
    else:
        _despike_input = vrAinput
    sl.append(
        "%s 3dDespike -nomask -overwrite -prefix %s %s "
        % (omp_threads_line_dspk_max, vrAinput, _despike_input)
    )
    vrAinput = "./%s_vrA%s" % (vrbase, osf)
if options.axialize:
    sl.append("3daxialize -overwrite -prefix %s %s" % (vrAinput, vrAinput))
    vrAinput = "./%s_vrA%s" % (vrbase, osf)
# Set eBbase
external_eBbase = False
if options.align_base != "":
    if options.align_base.isdigit():
        basevol = "%s[%s]" % (vrAinput, options.align_base)
    else:
        basevol = options.align_base
        external_eBbase = True
else:
    basevol = "%s[%s]" % (vrAinput, basebrik)
sl.append("3dcalc -a %s  -expr 'a' -prefix eBbase.nii.gz " % (basevol))
if external_eBbase:
    if oblique_mode:
        sl.append("3dWarp -wsinc5 -overwrite -deoblique eBbase.nii.gz eBbase.nii.gz")
    if options.axialize:
        sl.append("3daxialize -overwrite -prefix eBbase.nii.gz eBbase.nii.gz")
# Compute motion parameters
sl.append("volreg() {")
sl.append(
    "3dvolreg -overwrite -quintic  -prefix ./%s_vrA%s -base eBbase.nii.gz -dfile ./%s_vrA.1D -1Dmatrix_save ./%s_vrmat.aff12.1D %s"
    % (vrbase, osf, vrbase, PREFIX, vrAinput)
)
vrAinput = "./%s_vrA%s" % (vrbase, osf)
sl.append("1dcat './%s_vrA.1D[1..6]{%s..$}' > motion.1D " % (vrbase, basebrik))
sl.append("}")
sl.append(f"volreg {para_echo_indicator}")
e0dsin = PREFIX + DATASETS[0] + TRAILING

logcomment(
    "Preliminary preprocessing of functional datasets: despike, tshift, deoblique, and/or axialize",
    level=1,
)
# Do preliminary preproc for this run
if SHORTHAND_DSIN:
    DATASETS.sort()

para_proc_fun_names = []
for echo_ii in range(len(DATASETS)):
    # Determine dataset name
    echo = DATASETS[echo_ii]
    indata = getdsname(echo_ii)
    dsin = "e" + echo
    _para_proc_fun_name = f"prelim_loop_e{echo_ii}"
    sl.append(f"{_para_proc_fun_name} ()" + "{")
    logcomment(
        "Preliminary preprocessing dataset %s of TE=%sms to produce %s_ts+orig"
        % (indata, str(tes[echo_ii]), dsin)
    )
    # Pre-treat datasets: De-spike, RETROICOR in the future?
    intsname = "%s%s" % (dsprefix(indata), ISF)
    if not options.no_despike:
        intsname = "./%s_pt.nii.gz" % dsprefix(indata)
        sl.append(
            "%s 3dDespike -nomask -overwrite -prefix %s %s%s"
            % (omp_threads_line_dspk, intsname, dsprefix(indata), ISF)
        )
    # Time shift datasets
    if options.tpattern == "bids":
        # breakpoint()
        _src_dir = SRC_FOLDERS[echo_ii]
        json_file = os.path.join(_src_dir, f"{dsprefix(indata)}.json")
        with open(json_file) as bids_ifh:
            _bids_in = json.load(bids_ifh)
            if "SliceTiming" not in _bids_in.keys():
                tpat_opt = ""
            else:
                json_st_file = json_file + ".st.1D"
                st_extract_cmd = f"{sys.executable} -c \"import json; open('{json_st_file}', 'w').write(' '.join([str(tt) for tt in json.load(open('{_src_dir}/{json_file}'))['SliceTiming']]))\" "
                sl.append(st_extract_cmd)
                tpat_opt = f"-tpattern @{json_st_file}"
    elif options.tpattern != "":
        tpat_opt = " -tpattern %s " % options.tpattern
    else:
        tpat_opt = ""
    sl.append(
        "3dTshift -wsinc9 %s -prefix ./%s_ts+orig %s" % (tpat_opt, dsin, intsname)
    )
    # Force +orig label on dataset
    sl.append("3drefit -view orig %s_ts*HEAD" % (dsin))
    if oblique_mode and options.anat == "":
        sl.append(
            "3dWarp -wsinc5 -overwrite -deoblique -prefix ./%s_ts+orig ./%s_ts+orig"
            % (dsin, dsin)
        )
    # Axialize functional dataset
    if options.axialize:
        sl.append(
            "3daxialize  -overwrite -prefix ./%s_ts+orig ./%s_ts+orig" % (dsin, dsin)
        )
    if oblique_mode:
        sl.append("3drefit -deoblique -TR %s %s_ts+orig" % (options.TR, dsin))
    else:
        sl.append("3drefit -TR %s %s_ts+orig" % (options.TR, dsin))
    sl.append("}")
    para_proc_fun_names.append(_para_proc_fun_name)
for _fun_name in para_proc_fun_names:
    sl.append(f"{_fun_name} {para_echo_indicator}")
sl.append("wait")


# Compute T2*, S0, and OC volumes from raw data
logcomment(
    "Prepare T2* and S0 volumes for use in functional masking and (optionally) anatomical-functional coregistration (takes a little while).",
    level=1,
)
dss = DATASETS
dss.sort()
stackline = ""
para_proc_fun_names = []
for echo_ii in range(len(dss)):
    _para_proc_fun_name = f"vrA_loop_e{echo_ii}"
    sl.append(f"{_para_proc_fun_name} ()" + "{")
    echo = DATASETS[echo_ii]
    dsin = "e" + echo
    sl.append(
        "%s 3dAllineate -overwrite -final NN -NN -float -1Dmatrix_apply %s_vrmat.aff12.1D'{%i..%i}' -base eBbase.nii.gz -input %s_ts+orig'[%i..%i]' -prefix %s_vrA.nii.gz"
        % (
            omp_threads_line,
            PREFIX,
            int(basebrik),
            int(basebrik) + 20,
            dsin,
            int(basebrik),
            int(basebrik) + 20,
            dsin,
        )
    )
    stackline += " %s_vrA.nii.gz" % (dsin)
    sl.append("}")
    para_proc_fun_names.append(_para_proc_fun_name)
for _fun_name in para_proc_fun_names:
    sl.append(f"{_fun_name} {para_echo_indicator}")
sl.append("wait")


sl.append(f"t2smap ()" + "{")

sl.append("3dZcat -prefix basestack.nii.gz %s" % (stackline))
sl.append(
    "%s %s -d basestack.nii.gz -e %s"
    % (sys.executable, "/".join([meica_src_dir, "libmeica", "t2smap.py"]), options.tes)
)
sl.append("3dUnifize -prefix ./ocv_uni+orig ocv.nii")
sl.append(
    "3dSkullStrip -no_avoid_eyes -prefix ./ocv_ss.nii.gz -overwrite -input ocv_uni+orig"
)
sl.append(
    "3dcalc -overwrite -a t2svm.nii -b ocv_ss.nii.gz -expr 'a*ispositive(a)*step(b)' -prefix t2svm_ss.nii.gz"
)
sl.append(
    "3dcalc -overwrite -a s0v.nii -b ocv_ss.nii.gz -expr 'a*ispositive(a)*step(b)' -prefix s0v_ss.nii.gz"
)
if options.axialize:
    sl.append("3daxialize -overwrite -prefix t2svm_ss.nii.gz t2svm_ss.nii.gz")
    sl.append("3daxialize -overwrite -prefix ocv_ss.nii.gz ocv_ss.nii.gz")
    sl.append("3daxialize -overwrite -prefix s0v_ss.nii.gz s0v_ss.nii.gz")
sl.append("}")
sl.append(f"t2smap")

# Resume from here on
if options.resume:
    sl = []
    sl.append("export AFNI_DECONFLICT=OVERWRITE")

# Calculate affine anatomical warp if anatomical provided, then combine motion correction and coregistration parameters
if options.anat != "":
    # Copy in anatomical and make sure its in +orig space
    logcomment("Copy anatomical into ME-ICA directory and process warps", level=1)
    if oblique_mode:
        sl.append(
            "3dWarp -verb -card2oblique %s[0] -overwrite  -newgrid 1.000000 -prefix ./%s_ob.nii.gz %s/%s | grep  -A 4 '# mat44 Obliquity Transformation ::'  > %s_obla2e_mat.1D"
            % (vrinput, anatprefix, ANAT_DIR, nsmprage, PREFIX)  # type: ignore
        )
    sl.append("cp %s/%s* ." % (ANAT_DIR, nsmprage))  # type: ignore
    abmprage = nsmprage  # type: ignore
    refanat = nsmprage  # type: ignore
    if options.space:
        sl.append(f"anat_proc_st2()" + "{")
        sl.append("afnibinloc=`which 3dSkullStrip`")
        if "/" in options.space:
            sl.append('ll="%s"; templateloc=${ll%%/*}/' % options.space)
            options.space = options.space.split("/")[-1]
        else:
            sl.append("templateloc=${afnibinloc%/*}")
        atnsmprage = "%s_at.nii.gz" % (dsprefix(nsmprage))  # type: ignore
        if not dssuffix(nsmprage).__contains__("nii"):  # type: ignore
            sl.append(
                "3dcalc -float -a %s -expr 'a' -prefix %s.nii.gz"
                % (nsmprage, dsprefix(nsmprage))  # type: ignore
            )
        logcomment(
            "If can't find affine-warped anatomical, copy native anatomical here, compute warps (takes a while) and save in start dir. ; otherwise link in existing files"
        )
        sl.append(
            "if [ ! -e %s/%s ]; then thispwd=`pwd`; cd %s; @auto_tlrc -no_ss -init_xform AUTO_CENTER -base ${templateloc}/%s -input %s.nii.gz -suffix _at; cd $thispwd"
            % (ANAT_DIR, atnsmprage, ANAT_DIR, options.space, dsprefix(nsmprage))  # type: ignore
        )
        # sl.append("cp %s.nii %s" % (dsprefix(atnsmprage), ANAT_DIR))
        sl.append("gzip -f %s/%s.nii" % (ANAT_DIR, dsprefix(atnsmprage)))
        sl.append(
            "else if [ ! -e %s/%s ]; then ln -s %s/%s .; fi"
            % (ANAT_DIR, atnsmprage, ANAT_DIR, atnsmprage)
        )
        refanat = "%s/%s" % (ANAT_DIR, atnsmprage)
        sl.append("fi")
        sl.append(
            "3dcopy %s/%s.nii.gz %s"
            % (ANAT_DIR, dsprefix(atnsmprage), dsprefix(atnsmprage))
        )
        sl.append(
            "rm -f %s+orig.*; 3drefit -view orig %s+tlrc "
            % (dsprefix(atnsmprage), dsprefix(atnsmprage))
        )
        sl.append(
            "3dAutobox -overwrite -prefix ./abtemplate.nii.gz ${templateloc}/%s"
            % options.space
        )
        abmprage = "abtemplate.nii.gz"
        sl.append("}")
        sl.append(f"anat_proc_st2 {para_echo_indicator}")
        sl.append("wait")

        if options.qwarp:
            logcomment(
                "If can't find non-linearly warped anatomical, compute, save back; otherwise link"
            )
            nlatnsmprage = "%s_atnl.nii.gz" % (dsprefix(nsmprage))  # type: ignore
            sl.append("if [ ! -e %s/%s ]; then " % (ANAT_DIR, nlatnsmprage))
            logcomment(
                "Compute non-linear warp to standard space using 3dQwarp (get lunch, takes a while) "
            )
            sl.append(
                "3dUnifize -overwrite -GM -prefix ./%su.nii.gz %s/%s"
                % (dsprefix(atnsmprage), ANAT_DIR, atnsmprage)
            )
            sl.append(
                "3dQwarp -iwarp -overwrite -resample -useweight -blur 2 2 -duplo -workhard -base ${templateloc}/%s -prefix %s/%snl.nii.gz -source %s/%su.nii.gz"
                % (
                    options.space,
                    ANAT_DIR,
                    dsprefix(atnsmprage),
                    ANAT_DIR,
                    dsprefix(atnsmprage),
                )
            )
            sl.append("fi")
            sl.append(
                "if [ ! -e %s/%s ]; then ln -s %s/%s .; fi"
                % (ANAT_DIR, nlatnsmprage, ANAT_DIR, nlatnsmprage)
            )
            refanat = "%s/%snl.nii.gz" % (ANAT_DIR, dsprefix(atnsmprage))

    # Set anatomical reference for anatomical-functional co-registration
    if oblique_mode:
        alnsmprage = "./%s_ob.nii.gz" % (anatprefix)  # type: ignore
    else:
        alnsmprage = "%s/%s" % (ANAT_DIR, nsmprage)  # type: ignore
    if options.coreg_mode == "lp-t2s":
        ama_alnsmprage = alnsmprage
        if options.axialize:
            ama_alnsmprage = os.path.basename(alnsmprage)
            sl.append(
                "3daxialize -overwrite -prefix ./%s %s" % (ama_alnsmprage, alnsmprage)
            )
        t2salignpath = "libmeica/alignp_mepi_anat.py"
        sl.append(
            "%s %s -t t2svm_ss.nii.gz -a %s -p mepi %s"
            % (
                sys.executable,
                "/".join([meica_src_dir, t2salignpath]),
                ama_alnsmprage,
                options.align_args,
            )
        )
        sl.append(
            "cp alignp.mepi/mepi_al_mat.aff12.1D ./%s_al_mat.aff12.1D" % anatprefix  # type: ignore
        )
    elif options.coreg_mode == "aea":
        logcomment(
            "Using AFNI align_epi_anat.py to drive anatomical-functional coregistration "
        )
        sl.append("3dcopy %s ./ANAT_ns+orig " % alnsmprage)
        sl.append(
            "align_epi_anat.py -anat2epi -volreg off -tshift off -deoblique off -anat_has_skull no -save_script aea_anat_to_ocv.tcsh -anat ANAT_ns+orig -epi ocv_uni+orig -epi_base 0 %s"
            % (options.align_args)
        )
        sl.append("cp ANAT_ns_al_mat.aff12.1D %s_al_mat.aff12.1D" % (anatprefix))  # type: ignore
    if options.space:
        tlrc_opt = "%s/%s::WARP_DATA -I" % (ANAT_DIR, atnsmprage)  # type: ignore
        inv_tlrc_opt = "%s/%s::WARP_DATA" % (ANAT_DIR, atnsmprage)  # type: ignore
        sl.append(
            "cat_matvec -ONELINE %s > %s/%s_xns2at.aff12.1D"
            % (tlrc_opt, ANAT_DIR, anatprefix)  # type: ignore
        )
        sl.append(
            "cat_matvec -ONELINE %s > %s_xat2ns.aff12.1D" % (inv_tlrc_opt, anatprefix)  # type: ignore
        )
    else:
        tlrc_opt = ""
    if oblique_mode:
        oblique_opt = "%s_obla2e_mat.1D" % PREFIX
    else:
        oblique_opt = ""
    # pre-Mar 3, 2017, included tlrc affine warp in preprocessing. For new export flexiblity, will do tlrc_opt at export.
    # pre-Mar 3, 2017 version: sl.append("cat_matvec -ONELINE  %s %s %s_al_mat.aff12.1D -I > %s_wmat.aff12.1D" % (tlrc_opt,oblique_opt,anatprefix,prefix))
    sl.append(
        "cat_matvec -ONELINE  %s %s_al_mat.aff12.1D -I > %s_wmat.aff12.1D"
        % (oblique_opt, anatprefix, PREFIX)  # type: ignore
    )
    if options.anat:
        sl.append(
            "cat_matvec -ONELINE  %s %s_al_mat.aff12.1D -I  %s_vrmat.aff12.1D  > %s_vrwmat.aff12.1D"
            % (oblique_opt, anatprefix, PREFIX, PREFIX)  # type: ignore
        )
else:
    sl.append("cp %s_vrmat.aff12.1D %s_vrwmat.aff12.1D" % (PREFIX, PREFIX))


# Preprocess datasets
if SHORTHAND_DSIN:
    DATASETS.sort()

# Compute grand mean scaling factor
# sl.append(
#     "3dBrickStat -mask eBbase.nii.gz -percentile 50 1 50 e1_ts+orig[%i] > gms.1D"
#     % basebrik
# )
# sl.append("gms=`cat gms.1D`; gmsa=($gms); p50=${gmsa[1]}")

# Set resolution variables
sl.append(
    "voxsize=`ccalc .85*$(3dinfo -voxvol eBbase.nii.gz)**.33`"
)  # Set voxel size for decomp to slightly upsampled version of isotropic appx of native resolution so GRAPPA artifact is not at Nyquist
sl.append(
    'voxdims="`3dinfo -adi eBbase.nii.gz` `3dinfo -adj eBbase.nii.gz` `3dinfo -adk eBbase.nii.gz`"'
)
sl.append("echo $voxdims > voxdims.1D")
sl.append("echo $voxsize > voxsize.1D")

# Prepare eBvrmask
logcomment("Preparing functional masking for this ME-EPI run", level=1)
echo_ii = 0
echo = DATASETS[0]
indata = getdsname(0)
dsin = "e" + echo  # Note using same dsin as in time shifting
# abmprage = refanat = nsmprage   #Update as of Mar 3, 2017, to move to all native analysis
if options.anat:
    almaster = "-master %s" % nsmprage  # abmprage  # type: ignore
else:
    almaster = ""
# print 'almaster line is', almaster   #DEBUG
# print 'refanat line is', refanat   #DEBUG
sl.append("3dZeropad %s -prefix eBvrmask.nii.gz ocv_ss.nii.gz[0]" % (zeropad_opts))
# Create base mask
# if valid_qwarp_mode:
# 	if old_qwarp: nwarpstring = " -nwarp '%s/%s_WARP.nii.gz' -affter '%s_wmat.aff12.1D'" % (startdir,dsprefix(nlatnsmprage),prefix)
# 	else: nwarpstring = " -nwarp '%s/%s_WARP.nii.gz %s_wmat.aff12.1D' " % (startdir,dsprefix(nlatnsmprage),prefix)
if options.anat:
    sl.append(
        "3dAllineate -overwrite -final %s -%s -float -1Dmatrix_apply %s_wmat.aff12.1D -base %s -input eBvrmask.nii.gz -prefix ./eBvrmask.nii.gz %s %s"
        % ("NN", "NN", PREFIX, nsmprage, almaster, alfres)  # type: ignore
    )
    if options.t2salign or options.mask_mode != "func":
        sl.append(
            "3dAllineate -overwrite -final %s -%s -float -1Dmatrix_apply %s_wmat.aff12.1D -base eBvrmask.nii.gz -input t2svm_ss.nii.gz -prefix ./t2svm_ss_vr.nii.gz %s %s %s"
            % ("NN", "NN", PREFIX, almaster, alfres, para_echo_indicator)
        )
        sl.append(
            "3dAllineate -overwrite -final %s -%s -float -1Dmatrix_apply %s_wmat.aff12.1D -base eBvrmask.nii.gz -input ocv_uni+orig -prefix ./ocv_uni_vr.nii.gz %s %s %s"
            % ("NN", "NN", PREFIX, almaster, alfres, para_echo_indicator)
        )
        sl.append(
            "3dAllineate -overwrite -final %s -%s -float -1Dmatrix_apply %s_wmat.aff12.1D -base eBvrmask.nii.gz -input s0v_ss.nii.gz -prefix ./s0v_ss_vr.nii.gz %s %s %s"
            % ("NN", "NN", PREFIX, almaster, alfres, para_echo_indicator)
        )
        sl.append("wait")
# Fancy functional masking
if options.anat and options.mask_mode != "func":
    if options.space and options.mask_mode == "template":
        sl.append(
            "3dfractionize -overwrite -template eBvrmask.nii.gz -input abtemplate.nii.gz -prefix ./anatmask_epi.nii.gz -clip 1"
        )
        sl.append(
            "3dAllineate -overwrite -float -1Dmatrix_apply %s_xat2ns.aff12.1D -base eBvrmask.nii.gz -input anatmask_epi.nii.gz -prefix anatmask_epi.nii.gz -overwrite"
            % (anatprefix)  # type: ignore
        )
        logcomment(
            "Preparing functional mask using information from standard space template (takes a little while)"
        )
    if options.mask_mode == "anat":
        sl.append(
            "3dfractionize -template eBvrmask.nii.gz -input %s -prefix ./anatmask_epi.nii.gz -clip 0.5"
            % (nsmprage)  # type: ignore
        )
        logcomment(
            "Preparing functional mask using information from anatomical (takes a little while)"
        )
    sl.append(
        "3dBrickStat -mask eBvrmask.nii.gz -percentile 50 1 50 t2svm_ss_vr.nii.gz > t2s_med.1D"
    )
    sl.append(
        "3dBrickStat -mask eBvrmask.nii.gz -percentile 50 1 50 s0v_ss_vr.nii.gz > s0v_med.1D"
    )
    sl.append("t2sm=`cat t2s_med.1D`; t2sma=($t2sm); t2sm=${t2sma[1]}")
    sl.append("s0vm=`cat s0v_med.1D`; s0vma=($s0vm); s0vm=${s0vma[1]}")
    sl.append(
        '3dcalc -a ocv_uni_vr.nii.gz -b anatmask_epi.nii.gz -c t2svm_ss_vr.nii.gz -d s0v_ss_vr.nii.gz -expr "a-a*equals(equals(b,0)+isnegative(c-${t2sm})+ispositive(d-${s0vm}),3)" -overwrite -prefix ocv_uni_vr.nii.gz '
    )
    sl.append(
        "3dSkullStrip -no_avoid_eyes -overwrite -input ocv_uni_vr.nii.gz -prefix eBvrmask.nii.gz "
    )
    if options.fres:
        resstring = "-dxyz %s %s %s" % (
            options.fres,
            options.fres,
            options.fres,
        )
    else:
        resstring = "-dxyz ${voxsize} ${voxsize} ${voxsize}"
    sl.append(
        "3dresample -overwrite -master %s %s -input eBvrmask.nii.gz -prefix eBvrmask.nii.gz"
        % (nsmprage, resstring)  # type: ignore
    )

logcomment("Trim empty space off of mask dataset and/or resample")
sl.append("3dAutobox -overwrite -prefix eBvrmask%s eBvrmask%s" % (osf, osf))
resstring = "-dxyz ${voxsize} ${voxsize} ${voxsize}"
sl.append(
    "3dresample -overwrite -master eBvrmask.nii.gz %s -input eBvrmask.nii.gz -prefix eBvrmask.nii.gz"
    % (resstring)
)  # want this isotropic so spatial ops in select_model not confounded
sl.append(
    "3dcalc -float -a eBvrmask.nii.gz -expr 'notzero(a)' -overwrite -prefix eBvrmask.nii.gz"
)


logcomment("Extended preprocessing of functional datasets", level=1)
para_proc_fun_names = []
ica_input_list = []
for echo_ii in range(len(DATASETS)):
    # Determine dataset name
    _para_proc_fun_name = f"finalpp_loop_e{echo_ii}"
    sl.append(f"{_para_proc_fun_name} ()" + "{")
    echo = DATASETS[echo_ii]
    indata = getdsname(echo_ii)
    dsin = "e" + echo  # Note using same dsin as in time shifting
    # logcomment("Extended preprocessing dataset %s of TE=%sms to produce %s_in.nii.gz" % (indata,str(tes[echo_ii]),dsin),level=2 )
    logcomment(
        "Apply combined co-registration/motion correction parameter set to %s_ts+orig"
        % dsin
    )
    sl.append(
        "%s 3dAllineate -final %s -%s -float -1Dmatrix_apply %s_vrwmat.aff12.1D -base eBvrmask%s -input  %s_ts+orig -prefix ./%s_in%s"
        % (
            omp_threads_line,
            align_interp_final,
            align_interp,
            PREFIX,
            osf,
            dsin,
            dsin,
            osf,
        )
    )
    # sl.append("3dTstat -min -prefix ./%s_vr_min%s ./%s_vr%s" % (dsin, osf, dsin, osf))
    # sl.append(
    #     "3dcalc -a eBvrmask.nii.gz -b %s_vr_min%s -expr 'step(a)*step(b)' -overwrite -prefix eBvrmask.nii.gz "
    #     % (dsin, osf)
    # )
    if options.FWHM == "0mm":
        sl.append(
            "3dcalc -float -overwrite -a eBvrmask.nii.gz -b ./%s_in%s[%i..$] -expr 'step(a)*b' -prefix ./%s_in%s "
            % (dsin, osf, basebrik, dsin, osf)
        )
    else:
        sl.append(
            "3dBlurInMask -fwhm %s -mask eBvrmask%s -prefix ./%s_in%s ./%s_in%s[%i..$]"
            % (options.FWHM, osf, dsin, osf, dsin, osf, basebrik)
        )
    if options.detrend or options.highpass:
        sl.append("3dTstat -prefix ./%s_mean%s ./%s_in%s" % (dsin, osf, dsin, osf))
        if options.detrend:
            sl.append(
                "3dDetrend -polort %s -overwrite -prefix ./%s_in%s ./%s_in%s "
                % (options.detrend, dsin, osf, dsin, osf)
            )
        if options.highpass:
            sl.append(
                "3dBandpass -prefix ./%s_in%s %f 99 ./%s_in%s "
                % (dsin, osf, float(options.highpass), dsin, osf)
            )
        sl.append(
            "3dcalc -float -overwrite -a ./%s_in%s -b ./%s_mean%s -expr 'a+b' -prefix ./%s_in%s"
            % (dsin, osf, dsin, osf, dsin, osf)
        )
    ica_input_list.append(f"./{dsin}_in{osf}")
    # else:
    #     sl.append("mv %s_sm%s %s_in%s" % (dsin, osf, dsin, osf))
    # sl.append("3dTstat -stdev -prefix ./%s_std%s ./%s_in%s" % (dsin, osf, dsin, osf))
    sl.append("}")
    if not (options.test_proc or options.keep_int):
        sl.append("rm -f %s_pt.nii.gz %s_vr%s %s_sm%s" % (dsin, dsin, osf, dsin, osf))
    para_proc_fun_names.append(_para_proc_fun_name)
for _fun_name in para_proc_fun_names:
    sl.append(f"{_fun_name} {para_echo_indicator}")
sl.append("wait")

ica_input = " ".join(ica_input_list)

# Spatial concatenation of datasets, this needs to get removed in future versions based on the argparse feature.
# ica_input = "zcat_ffd.nii.gz"
# ica_mask = "zcat_mask.nii.gz"
# zcatstring = ""
# for echo in ica_datasets:
#     dsin = "e" + echo
#     zcatstring = "%s ./%s_in%s" % (zcatstring, dsin, osf)
# sl.append("3dZcat -overwrite -prefix %s  %s" % (ica_input, zcatstring))
# sl.append(
#     "3dcalc -float -overwrite -a %s[0] -expr 'notzero(a)' -prefix %s"
#     % (ica_input, ica_mask)
# )

if options.pp_only:
    tedflag = "#"
else:
    tedflag = ""

if options.resume:
    sl.append("rm -f TED/pcastate.pklbz")
if options.tedica_only:
    sl = []

strict_setting = ""
if options.strict:
    strict_setting = "--strict"

ted_spms_setting = ""
if options.ted_spms:
    ted_spms_setting = "--fout"

# if os.path.exists("%s/libmeica" % (meicadir)):
tedanapath = "tedana.py"
# else:
#     tedanapath = "lib_tedana.py"

ted_seed_setting = ""
if options.ted_seed:
    ted_seed_setting = f"--seed={options.ted_seed}"

logcomment("Perform TE-dependence analysis (takes a good while)", level=1)
interactive_flag = ""  # deprecated

T2arg = "T2brain=80"
tedcall = (
    "STARTDIR=%s/meica.%s %s%s %s %s -e %s  -d %s --sourceTEs=%s --kdaw=%s --rdaw=%s --initcost=%s --finalcost=%s --post_gscontrol --pre_gscontrol --conv=2.5e-5 %s %s %s %s %s"
    % (
        STARTDIR,
        SETNAME,
        tedflag,
        sys.executable,
        interactive_flag,
        "/".join([meica_src_dir, tedanapath]),
        options.tes,
        ica_input,
        options.sourceTEs,
        options.daw,
        "1",
        options.initcost,
        options.finalcost,
        strict_setting,
        ted_spms_setting,
        ted_seed_setting,
        options.ted_args,
        T2arg,
    )
)

tedcall_debug_cs = tedcall + " --mix=meica_mix.1D DEBUG"
sl.append(f"echo {tedcall_debug_cs} > debug_component_selection.sh")

sl.append(tedcall)

if outprefix == "":
    outprefix = SETNAME

if options.select_only:
    sl = []
    # sl.append(
    #     "voxsize=`ccalc .85*$(3dinfo -voxvol eBbase.nii.gz)**.33`"
    # )  # Set voxel size for decomp to slightly upsampled version of isotropic appx of native resolution so GRAPPA artifact is not at Nyquist
    # sl.append(
    #     'voxdims="`3dinfo -adi eBbase.nii.gz` `3dinfo -adj eBbase.nii.gz` `3dinfo -adk eBbase.nii.gz`"'
    # )
    # sl.append("echo $voxdims > voxdims.1D")
    # sl.append("echo $voxsize > voxsize.1D")
    sl.append(
        "STARTDIR=%s/%s %s %s %s -e %s  -d %s --sourceTEs=%s --mix=meica_mix.1D --post_gscontrol --pre_gscontrol %s %s %s"
        % (
            STARTDIR,
            SETNAME,
            sys.executable,
            interactive_flag,
            "/".join([meica_src_dir, tedanapath]),
            options.tes,
            ica_input,
            options.sourceTEs,
            options.ted_args,
            VENDOR_OPT,
            VENDOR_ARGS,
        )
    )

mask_dict = {}  # Need this here


def export_glob(
    inglob,
    outfileprefix,
    comment="Created by %s" % runcmd,
    interp="wsinc5",
    disable=False,
    plaintext=False,
    export_mask="./export_mask.nii.gz",
    globsuffix=None,
):
    meicasuffix = None
    exportmode = None

    # Compile list of files to transform in one call
    sl.append('FILE_LIST=""')
    sl.append(f"for FILE in `ls {inglob}`; do")
    sl.append('    FILE_LIST="${FILE_LIST} ${FILE}";')
    sl.append("done")

    # Transform in one step if Nwarp
    if valid_qwarp_mode:
        exportmode = "WARPONLY"
        export_result(
            "${FILE_LIST}",
            outfileprefix=None,
            comment=comment,
            interp=interp,
            disable=disable,
            plaintext=plaintext,
            export_mask=export_mask,
            exportmode=exportmode,
            globsuffix=globsuffix,
        )
        exportmode = "POSTONLY"
        meicasuffix = "_nlw"

    # If file is 'foo/bar/baz.nii.gz'
    #   and globsuffix is _ted
    # FILE_DIR is foo/bar
    # FILE_BASE is baz.nii.gz
    # FILE_PREFIX is baz
    # FILE_TYPE is .nii
    # INFILE is baz_ted.nii
    #
    sl.append(f"for FILE in $FILE_LIST; do")
    sl.append("    FILE_DIR=${FILE%/*}")
    sl.append("    FILE_BASE=${FILE#*/}")
    sl.append("    FILE_PREFIX=${FILE_BASE%%.*}")
    sl.append("    FILE_SUFFIX=%s" % globsuffix)
    sl.append("    FILE_TYPE=.nii")
    if valid_qwarp_mode:
        sl.append("    INFILE=${FILE_BASE%%.*}${FILE_SUFFIX}${FILE_TYPE}")
        sl.append(
            "mv ${FILE_PREFIX}%s${FILE_TYPE} %s${FILE_PREFIX}%s%s${FILE_TYPE}"
            % (globsuffix, outfileprefix, globsuffix, meicasuffix)
        )
    else:
        exportmode = "WARPANDPOST"
        sl.append("    INFILE=$FILE_BASE")
    export_result(
        "${FILE_DIR}/${INFILE}",
        outfileprefix="%s${FILE_PREFIX}%s" % (outfileprefix, globsuffix),
        comment=comment,
        interp=interp,
        disable=disable,
        plaintext=plaintext,
        export_mask=export_mask,
        exportmode=exportmode,  # type: ignore
        globsuffix="",
    )
    sl.append("done")


def export_result(
    infile,
    outfileprefix,
    comment="Created by %s" % runcmd,
    interp="wsinc5",
    disable=False,
    plaintext=False,
    export_mask="./export_mask.nii.gz",
    exportmode="WARPANDPOST",
    globsuffix=None,
):
    logcomment(f"Export mode is {exportmode}", 4)

    esl = []
    tedflag = ""
    if disable == True:
        tedflag = "#"
    native_export = options.native
    if plaintext:
        sl.append("cp %s %s/%s;" % (infile, STARTDIR, outfileprefix))
        return
    # Set output resolution parameters, either specified fres or original voxel dimensions
    if options.fres:
        resstring = "-dxyz %s %s %s" % (options.fres, options.fres, options.fres)
    else:
        resstring = "-dxyz ${voxdims}"
    to_export = []
    if options.space:
        export_master = "-master %s" % abmprage  #

    # If Qwarp, do Nwarpapply
    if valid_qwarp_mode:
        warp_code = "nlw"

        suffixopt = ""
        prefixopt = ""
        if globsuffix is not None:
            suffixopt = f"-suffix {globsuffix}"
        if outfileprefix is not None:
            prefixopt = f"-prefix {outfileprefix}_{warp_code}.nii"

        this_nwarpstring = " -nwarp %s/%s_xns2at.aff12.1D '%s/%s_WARP.nii.gz' " % (
            ANAT_DIR,
            anatprefix,
            STARTDIR,
            dsprefix(nlatnsmprage),
        )
        esl.append(
            "%s3dNwarpApply -overwrite %s %s %s -source %s -interp %s %s %s;"
            % (
                tedflag,
                this_nwarpstring,
                export_master,  # type: ignore
                qwfres,
                infile,
                interp,
                suffixopt,
                prefixopt,
            )
        )
        if not warp_code in list(mask_dict.keys()):
            prefixopt = f"-prefix {warp_code}_export_mask.nii"
            esl.append(
                "%s3dNwarpApply -overwrite %s %s %s -source %s -interp %s %s;"
                % (
                    tedflag,
                    this_nwarpstring,
                    export_master,  # type: ignore
                    qwfres,
                    "export_mask.nii.gz",
                    interp,
                    prefixopt,
                )
            )
            esl.append(
                "%snifti_tool -mod_hdr -mod_field sform_code 2 -mod_field qform_code 2 -infiles %s_export_mask.nii -overwrite;"
                % (tedflag, warp_code)
            )
            mask_dict[warp_code] = "%s_export_mask.nii" % warp_code
        to_export.append(
            ("%s_%s" % (outfileprefix, warp_code), "%s_export_mask.nii" % warp_code)
        )
    # If there's a template space, allineate result to that space
    elif options.space:
        warp_code = "afw"
        esl.append(
            "%s3dAllineate -overwrite -final %s -%s -float -1Dmatrix_apply %s/%s_xns2at.aff12.1D -input %s -prefix ./%s_afw.nii %s %s;"
            % (
                tedflag,
                interp,
                align_interp,
                ANAT_DIR,
                anatprefix,
                infile,
                outfileprefix,
                export_master,  # type: ignore
                alfres,
            )
        )
        if not warp_code in list(mask_dict.keys()):
            esl.append(
                "%s3dAllineate -overwrite -final %s -%s -float -1Dmatrix_apply %s/%s_xns2at.aff12.1D -input %s -prefix ./%s_export_mask.nii %s %s;"
                % (
                    tedflag,
                    interp,
                    align_interp,
                    ANAT_DIR,
                    anatprefix,
                    "export_mask.nii.gz",
                    warp_code,
                    export_master,  # type: ignore
                    alfres,
                )
            )
            esl.append(
                "%snifti_tool -mod_hdr -mod_field sform_code 2 -mod_field qform_code 2 -infiles %s_export_mask.nii -overwrite;"
                % (tedflag, warp_code)
            )
            mask_dict[warp_code] = "%s_export_mask.nii" % warp_code
        to_export.append(
            ("%s_%s" % (outfileprefix, warp_code), "%s_export_mask.nii" % warp_code)
        )
    # Otherwise resample
    else:
        native_export = True
    if native_export:
        if options.anat:
            native_suffix = "nat"
            export_master = "-master %s" % nsmprage
        else:
            native_suffix = "epi"
            export_master = ""
        esl.append(
            "%s3dresample -rmode Li -overwrite %s %s -input %s -prefix %s_%s.nii;"
            % (tedflag, export_master, resstring, infile, outfileprefix, native_suffix)
        )
        warp_code = native_suffix
        if not warp_code in list(mask_dict.keys()):
            # sl.append('3dresample -overwrite -prefix %s_export_mask.nii -rmode Li -master %s_%s.nii -input export_mask.nii.gz' % (warp_code,outfileprefix,warp_code))
            esl.append(
                "%s3dresample -rmode Li -overwrite %s %s -input %s -prefix %s_export_mask.nii;"
                % (tedflag, export_master, resstring, "export_mask.nii.gz", warp_code)
            )
            mask_dict[warp_code] = "%s_export_mask.nii" % warp_code
        to_export.append(
            ("%s_%s" % (outfileprefix, warp_code), "%s_export_mask.nii" % warp_code)
        )
    if exportmode == "WARPANDPOST" or exportmode == "WARPONLY":
        sl.extend(esl)
    if exportmode == "WARPANDPOST" or exportmode == "POSTONLY":
        for P, P_mask in to_export:
            export_postproc(P, P_mask, tedflag, comment, STARTDIR)


def export_postproc(P, P_mask, tedflag, comment, startdir):
    # sl.append("%s3dNotes -h '%s' %s.nii;" % (tedflag, comment, P))
    if options.anat != "" and options.space != False and "_nat" not in P:
        sl.append(
            "%snifti_tool -mod_hdr -mod_field sform_code 2 -mod_field qform_code 2 -infiles %s.nii -overwrite;"
            % (tedflag, P)
        )
    sl.append(
        "3dcalc -overwrite -a %s -b %s.nii -expr 'ispositive(a-.5)*b' -prefix %s.nii.gz ; mv %s.nii.gz %s;"
        % (P_mask, P, P, P, startdir)
    )


if options.export_only:
    sl = []

sl.append(
    'voxdims="`3dinfo -adi eBbase.nii.gz` `3dinfo -adj eBbase.nii.gz` `3dinfo -adk eBbase.nii.gz`"'
)
sl.append("echo $voxdims > voxdims.1D")
# Make the export mask
sl.append(
    "3dcalc -float -a TED/ts_OC.nii[0] -overwrite -expr 'notzero(a)' -prefix ./export_mask.nii.gz"
)
logcomment("Copying results to start directory", level=1)
export_result(
    "TED/ts_OC.nii",
    "%s_tsoc" % (outprefix),
    "T2* weighted average of ME time series, produced by ME-ICA %s" % __version__,
    disable=options.pp_only,
)
export_result(
    "TED/dn_ts_OC.nii",
    "%s_medn" % (outprefix),
    "Denoised timeseries (including thermal noise), produced by ME-ICA %s"
    % __version__,
    disable=options.pp_only,
)
export_result(
    "TED/dn_ts_OC_T1c.nii",
    "%s_T1c_medn" % (outprefix),
    "Denoised timeseries with T1 equilibration correction (including thermal noise), produced by ME-ICA %s"
    % __version__,
    disable=options.pp_only,
)
export_result(
    "TED/hik_ts_OC_T1c.nii",
    "%s_hikts" % (outprefix),
    "Denoised timeseries with T1 equilibration correction (no thermal noise), produced by ME-ICA %s"
    % __version__,
    disable=options.pp_only,
)
export_result(
    "TED/hik_ts_OC_T1c.dr2s.nii",
    "%s_dr2s" % (outprefix),
    "Denoised timeseries with T1 equilibration correction (no thermal noise), produced by ME-ICA %s"
    % __version__,
    disable=options.pp_only,
)
export_result(
    "TED/betas_hik_OC.nii",
    "%s_mefc" % (outprefix),
    "Denoised ICA coeff. set for ME-ICR seed-based FC analysis, produced by ME-ICA %s"
    % __version__,
    disable=options.pp_only,
)
export_result(
    "TED/betas_OC.nii",
    "%s_mefl" % (outprefix),
    "Full ICA coeff. set for component assessment, produced by ME-ICA %s" % __version__,
    disable=options.pp_only,
)
export_result(
    "TED/feats_OC2.nii",
    "%s_mefcz" % (outprefix),
    "Z-normalized spatial component maps, produced by ME-ICA %s" % __version__,
    disable=options.pp_only,
)
export_result(
    "TED/t2svG.nii",
    "%s_t2s" % (outprefix),
    "Proprietary load metrics produced by ME-ICA %s" % __version__,
    disable=options.pp_only,
)
export_result("TED/comp_table.txt", "%s_ctab.txt" % (outprefix), plaintext=True)
export_result("TED/meica_mix.1D", "%s_mmix.1D" % (outprefix), plaintext=True)
if options.ted_spms:
    # Output TE-dependence maps
    export_glob(
        "TED/cc*.nii",
        outfileprefix=f"{outprefix}",
        globsuffix="_ted",
        comment="TE-dependence map, produced by ME-ICA %s" % __version__,
        disable=options.pp_only,
    )

# VENDOR OUTPUTS
if VENDOR_OPT:
    v_inputs, v_outputs, v_descs = vendor_exports(outprefix)
    for _ei, _input in enumerate(v_inputs):
        export_result(
            _input,
            v_outputs[_ei],
            v_descs[_ei],
            disable=options.pp_only,
        )


if initial_run and not options.no_trim_work:
    e_csl = ",".join([str(ii) for ii in list(range(1, len(DATASETS) + 1))])
    sl.append(f"mv {workdir_full_path}/e[{e_csl}]_in{osf} {workdir_full_path}/TED/")
    sl.append(f"rm {workdir_full_path}/e[{e_csl}]*nii*")
    sl.append(f"rm {workdir_full_path}/e[{e_csl}]*BRIK*")
    sl.append(f"rm {workdir_full_path}/e[{e_csl}]*HEAD*")
    sl.append(f"rm {workdir_full_path}/*_pt*")
    sl.append(f"rm {workdir_full_path}/*_vrA*")
    sl.append(f"rm {workdir_full_path}/*.nii")
    sl.append(f"rm -rf {workdir_full_path}/memmaps")
    sl.append(f"mv {workdir_full_path}/TED/e[{e_csl}]_in{osf} {workdir_full_path}/")
    sl.append(f"cp {workdir_full_path}/TED/meica_mix.1D {workdir_full_path}/")
    sl.append(f"cp {workdir_full_path}/TED/encoding_fspace*.1D {workdir_full_path}/")
    sl.append(f"cp {workdir_full_path}/TED/comp_table.txt {workdir_full_path}/")
    sl.append(f"cp {workdir_full_path}/TED/betas_OC.nii {workdir_full_path}/")
    sl.append(f"rm -rf {workdir_full_path}/TED/*")
    sl.append(f"mv {workdir_full_path}/comp_table.txt {workdir_full_path}/TED/")
    sl.append(f"mv {workdir_full_path}/meica_mix.1D {workdir_full_path}/TED/")
    sl.append(f"mv {workdir_full_path}/encoding_fspace*.1D {workdir_full_path}/TED/")
    sl.append(f"mv {workdir_full_path}/betas_OC.nii {workdir_full_path}/TED/")

if linked_workdir and initial_run:
    sl.append(f"rm {finaldir_full_path}")
    # replace with minimal file copy later
    sl.append(f"cp -rf {workdir_full_path} {finaldir_full_path}")
    sl.append(f"rm -rf {workdir_full_path}")

# Write the preproc script and execute it
ofh = open(f"{STARTDIR}/_meica_%s.sh" % SETNAME, "w")
print("++ Writing script file: _meica_%s.sh" % (SETNAME))
ofh.write("\n".join(headsl + sl) + "\n")
ofh.close()
if not options.script_only:
    # should add suffixes for only so we don't overwrite origina _meica...sh file
    print("++ Executing script file: _meica_%s.sh" % (SETNAME))
    system("bash _meica_%s.sh" % SETNAME)
