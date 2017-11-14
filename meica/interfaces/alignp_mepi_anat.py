#!/usr/bin/env python

import re
import sys
import os.path
import subprocess
from re import split as resplit
from os import system, getcwd, mkdir, chdir, popen


def dsprefix(idn):
    def prefix(datasetname):
        return datasetname.split('+')[0]
    if len(idn.split('.'))!=0:
        if idn.split('.')[-1]=='HEAD' or idn.split('.')[-1]=='BRIK' or idn.split('.')[-2:]==['BRIK','gz']:
            return prefix(idn)
        elif idn.split('.')[-1]=='nii' and not idn.split('.')[-1]=='nii.gz':
            return '.'.join(idn.split('.')[:-1])
        elif idn.split('.')[-2:]==['nii','gz']:
            return '.'.join(idn.split('.')[:-2])
        else:
            return prefix(idn)
    else:
        return prefix(idn)


def dssuffix(idna):
    suffix = idna.split(dsprefix(idna))[-1]
    spl_suffix=suffix.split('.')
    if len(spl_suffix[0])!=0 and spl_suffix[0][0] == '+': return spl_suffix[0]
    else: return suffix


def niibrik(nifti): return dsprefix(nifti)+'+orig'


def import_datasets(dsets):
    """
    Base is functional
    Dset is anatomical
    """
    outnames = []
    for dset in dsets:
        if dset!='' and dset!=None:
            indir=''
            if dset[0]!='/': indir = startdir
            sl.append("3dcopy -overwrite %s/%s %s/%s/%s" % (indir,dset,startdir,walignp_dirname,dsprefix(os.path.basename(dset))))
            sl.append("3drefit -view orig `ls %s/%s/%s+*.HEAD`" % (startdir,walignp_dirname,dsprefix(os.path.basename(dset))))
            impdset = '%s+orig' % dsprefix(os.path.basename(dset))
            outnames.append(niibrik(impdset))
        else:
            outnames.append('')
    return outnames


def align_centers(base,dset):
    sl.append("@Align_Centers -base %s+orig -dset %s+orig" % (dsprefix(base),dsprefix(dset)))
    volname = dsprefix(dset) + '_shft+orig'
    matname = dsprefix(dset) + '_shft.1D'
    return volname,matname


def makeflat(dset):
    flatname = "%s_flat" % (dsprefix(dset))
    sl.append("3dcalc -a %s -expr 'step(a)' -prefix %s" % (dset,flatname) )
    return niibrik(flatname)


def graywt(t2sname,s0name):
    basevol='align_base.nii.gz'
    weightvol='align_weight.nii.gz'
    if s0name!='':
        sl.append("3dcalc -overwrite -a %s -b %s -expr 'a*ispositive(a)*step(b)' -prefix t2svm_ss.nii.gz" % (t2sname,s0name) )
        t2sname='t2svm_ss.nii.gz'
    sl.append("3dBrickStat -mask %s -percentile 50 1 50 %s > graywt_thr.1D" %  (t2sname,t2sname) )
    sl.append("maxthr=`cat graywt_thr.1D`; maxthra=($maxthr); vmax=${maxthra[1]}" )
    sl.append("3dcalc -overwrite -a %s -expr \"a*isnegative(a-2*${vmax})\" -prefix t2svm_thr.nii.gz" % (t2sname) )
    sl.append("3dUnifize -overwrite -prefix align_base.nii.gz t2svm_thr.nii.gz")
    sl.append("3dSeg -prefix Segsy.t2svm -anat %s -mask %s" % (basevol,basevol))
    sl.append("3dcalc -overwrite -a 'Segsy.t2svm/Posterior+orig[2]' -prefix %s -expr 'a' " % weightvol)
    return basevol,weightvol


def allineate(sourcevol, weight, targetvol, prefix, maxrot, maxshf,maxscl,do_cmass):
    """
    source should be anatomical
    target should be T2*
    """
    outvol_prefix = "%s_al"  % prefix
    outmat_prefix = "%s_al_mat"  % prefix
    cmass_opt = ''
    autocmass=False
    if do_cmass or options.autocmass: cmass_opt = '-cmass'
    elif not do_cmass and not options.autocmass: cmass_opt=''
    align_opts = "-lpc -weight %s -maxshf %s -maxrot %s -maxscl %s %s" % (weightvol, maxrot, maxshf,maxscl,cmass_opt)
    sl.append("3dAllineate -overwrite -weight_frac 1.0 -VERB -warp aff -source_automask+2 -master SOURCE -source %s -base %s -prefix ./%s -1Dmatrix_save ./%s %s " \
        % (sourcevol,targetvol,outvol_prefix,outmat_prefix,align_opts))
    if options.autocmass:
        cmass_opt=''
        align_opts = "-lpc -weight %s -maxshf %s -maxrot %s -maxscl %s %s" % (weightvol, maxrot, maxshf,maxscl,cmass_opt)
    sl.append("3dAllineate -overwrite -weight_frac 1.0 -VERB -warp aff -source_automask+2 -master SOURCE -source %s -base %s -prefix ./ncm_%s -1Dmatrix_save ./ncm_%s %s " \
        % (sourcevol,targetvol,outvol_prefix,outmat_prefix,align_opts))
    outvol_name = niibrik(outvol_prefix)
    outmat_name = outmat_prefix + 'blah'

    return outvol_name, outmat_name


def cmselect(base,allin):
    if options.autocmass:
        sl.append("3dresample -overwrite -master %s -prefix rs_%s -inset %s -rmode NN" % (base,allin,allin) )
        sl.append("3dresample -overwrite -master %s -prefix rs_ncm_%s -inset ncm_%s -rmode NN" % (base,allin,allin) )
        sl.append("echo `3ddot -dodice %s rs_ncm_%s` > ncm_olap.txt" % (base,allin))
        sl.append("echo `3ddot -dodice %s rs_%s` > cm_olap.txt" % (base,allin))
        sl.append("ncm_olap=`cat ncm_olap.txt`; cm_olap=`cat cm_olap.txt`")
        sl.append("""mv_dec=`echo "$ncm_olap>$cm_olap" | bc`""")
        sl.append("""if [ "$mv_dec" == "1" ] ; then 3drename -overwrite ncm_%s %s ; mv ncm_%s_mat.aff12.1D %s_mat.aff12.1D; fi """ % (allin,allin,dsprefix(allin),dsprefix(allin)))


def runproc(script_prefix, scrlines):
    script_name = "_walignp_mepi_anat_%s.sh" % script_prefix
    ofh = open(script_name,'w')
    ofh.write("\n".join(scrlines))
    ofh.close()
    os.system("bash %s" % script_name)
