#!/usr/bin/env python

import os
import os.path as op
import numpy as np
import meica
import pytest

# change into the resources directory
os.chdir(op.join(op.dirname(op.abspath(__file__)), 'resources'))

# set file names for testing
f1 = 'sub-001_task-rest_run-01_echo-1_bold.nii.gz'
f2 = 'sub-001_task-rest_run-01_echo-2_bold.nii.gz'
f3 = 'sub-001_task-rest_run-01_echo-3_bold.nii.gz'


def test_fparse():
    # identify prefix, trailing file part, and file type for nii file type
    assert ('sub-001_task-rest_run-01_echo-',
            '_bold', 
            '.nii.gz') == meica.fparse([f1, f2, f3])

    # identify prefix, trailing file part, and file type for AFNI file type
    assert ('sub.001_e0',
            '_localizer',
            '+tlrc')   == meica.fparse(['sub.001_e01_localizer+tlrc.BRIK.gz', 
                                        'sub.001_e02_localizer+tlrc.BRIK.gz'])


def test_format_inset():
    # list dataset specification, string TE specification
    assert (f2,'sub-001_task-rest_run-01_echo-123_bold') == meica.format_inset([f1,f2,f3],
                                                                                '12.2,24.6,30')

    # string dataset specification, list TE specification
    assert (f2,'sub-001_task-rest_run-01_echo-123_bold') == meica.format_inset(','.join([f1, f2 ,f3]),
                                                                               [12.2, 24.6, 30])

    # shorthand dataset specification, string TE specification
    assert (f2, 'sub-001_task-rest_run-01_echo-123_bold') == meica.format_inset('sub-001_task-rest_run-01_echo-[1,2,3]_bold.nii.gz',
                                                                                '12.2,24.6,30')

    # catch differing number of datasets and echoes
    with pytest.raises(AssertionError):
        meica.format_inset(','.join([f1, f2 ,f3]),
                           [12.2, 24.6])

    # catch incorrect echo input
    with pytest.raises(TypeError):
        meica.format_tes(12.2)

    # catch incorrect dataset input
    with pytest.raises(TypeError):
        meica.format_dset(3)


    # ('sub_001.e02+tlrc.HEAD',
    #  'sub_001.e0123') == format_inset('sub_001.e0[1,2,3]+tlrc.BRIK.gz',
    #                                   [12.2, 24.6, 30])

    # ('sub_001.e02+tlrc.BRIK.gz',
    #  'sub_001.e0123') == format_inset('sub_001.e01+tlrc.BRIK.gz, ' +
    #                                   'sub_001.e02+tlrc.BRIK.gz, ' +
    #                                   'sub_001.e03+tlrc.BRIK.gz',
    #                                   [12.2, 24.6, 30])


def test_find_CM():
    assert np.allclose([1.1718978881835938,
                        -41.01300048828125,
                        -46.293296813964844], 
                        meica.find_CM('sub-001_T1w.nii.gz'))

# def test_gen_script():
#     os.environ['DYLD_FALLBACK_LIBRARY_PATH'] = '~/abin'
#     resdir = op.join(op.dirname(__file__),'resources')
#     fname = op.join(resdir,
#                     '_meica_sub-001_task-rest_echo-123_run-01_meepi.sh')
#     sel_opts = ['-d', op.join(resdir,
#                               'sub-001_task-rest_echo-[1,2,3]_run-01_meepi.nii.gz'),
#                 '-e', '14.5,38.5,62.5',
#                 '-b', '4v',
#                 '-a', op.join(resdir,'sub-001_T1w.nii.gz'),
#                 '--fres=2', '--MNI', '--qwarp']

#     opts = meica.get_options(_debug=sel_opts)
#     with open(fname, 'r') as file:
#         script = file.read()
#     script_list = meica.gen_script(opts)
#     "\n".join(script_list) == script
