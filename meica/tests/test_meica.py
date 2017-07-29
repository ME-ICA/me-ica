#!/usr/bin/env python

import os
import os.path as op
from meica.meica import find_CM, format_inset, fparse, gen_script, get_options

def test_fparse():
    ('sub.001', '.nii.gz') == fparse('sub.001.nii.gz')
    ('sub.001', '+tlrc')   == fparse('sub.001+tlrc.BRIK.gz')


def test_format_inset():
    ('sub_001.e02.nii.gz',
     'sub_001.e0123') == format_inset('sub_001.e0[1,2,3].nii.gz')

    ('sub_001.e02.nii.gz',
     'sub_001.e0123') == format_inset('sub_001.e01.nii.gz, ' +
                                      'sub_001.e02.nii.gz, ' +
                                      'sub_001.e03.nii.gz',
                                      [12.2, 24.6, 30])

    ('sub_001.e02+tlrc.HEAD',
     'sub_001.e0123') == format_inset('sub_001.e0[1,2,3]+tlrc.BRIK.gz',
                                      [12.2, 24.6, 30])

    ('sub_001.e02+tlrc.BRIK.gz',
     'sub_001.e0123') == format_inset('sub_001.e01+tlrc.BRIK.gz, ' +
                                      'sub_001.e02+tlrc.BRIK.gz, ' +
                                      'sub_001.e03+tlrc.BRIK.gz',
                                      [12.2, 24.6, 30])


def test_find_CM():
    resdir = op.join(op.dirname(__file__),'resources')

    [1.1718978881835938,
     -41.01300048828125,
     -46.293296813964844] == find_CM(op.join(resdir,
                                             '/sub-001_T1w.nii.gz'))

def test_gen_script():
    os.environ['DYLD_FALLBACK_LIBRARY_PATH'] = '~/abin'
    resdir = op.join(op.dirname(__file__),'resources')
    fname = op.join(resdir,
                    '_meica_sub-001_task-rest_echo-123_run-01_meepi.sh')
    sel_opts = ['-d', op.join(resdir,
                              'sub-001_task-rest_echo-[1,2,3]_run-01_meepi.nii.gz'),
                '-e', '14.5,38.5,62.5',
                '-b', '4v',
                '-a', op.join(resdir,'sub-001_T1w.nii.gz'),
                '--fres=2', '--MNI', '--qwarp']

    opts = get_options(_debug=sel_opts)
    with open(fname, 'r') as file:
        script = file.read()
    script_list = gen_script(opts)
    "\n".join(script_list) == script
