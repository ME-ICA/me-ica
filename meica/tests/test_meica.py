#!/usr/bin/env python

import os
import os.path as op
import numpy as np
import meica
from nipype.interfaces.traits_extension import TraitError
import pytest

try:
    FileNotFoundError
except NameError: # working in python 2
    FileNotFoundError = IOError

# change into the resources directory
os.chdir(op.join(op.dirname(op.abspath(__file__)), 'resources'))

# set file names for testing
nii1  = 'sub-001_task-rest_run-01_echo-1_bold.nii.gz'
nii2  = 'sub-001_task-rest_run-01_echo-2_bold.nii.gz'
nii3  = 'sub-001_task-rest_run-01_echo-3_bold.nii.gz'

afni1 = 'sub-001_task-rest_run-02_echo-1_bold+orig.BRIK.gz'
afni2 = 'sub-001_task-rest_run-02_echo-2_bold+orig.BRIK.gz'
afni3 = 'sub-001_task-rest_run-02_echo-3_bold+orig.BRIK.gz'


def test_afni_fname_parse():
    """
    """

    # identify prefix, trailing file part, and file type for AFNI file type
    assert ('sub-001_task-rest_run-02_echo-',
            '_bold',
            '+orig')   == meica.afni_fname_parse(afni1, echo_ind=30)


def test_nii_fname_parse():
    """
    """

    # identify prefix, trailing file part, and file type for nii file type
    assert ('sub-001_task-rest_run-01_echo-',
            '_bold',
            '.nii.gz') == meica.nii_fname_parse(nii1, echo_ind=30)


def test_fname_parse():
    """
    """

    # identify prefix, trailing file part, and file type for AFNI file type
    assert ('sub-001_task-rest_run-02_echo-',
            '_bold',
            '+orig')   == meica.fname_parse([afni1, afni2])

    # identify prefix, trailing file part, and file type for nii file type
    assert ('sub-001_task-rest_run-01_echo-',
            '_bold',
            '.nii.gz') == meica.fname_parse([nii1, nii2, nii3])

    # identify prefix, trailing file part, and file type for anatomical
    assert ('sub-001_T1w',
            '',
            '.nii.gz') == meica.fname_parse('sub-001_T1w.nii.gz')

    # catch incorrect file list input
    with pytest.raises(TypeError) as err:
        meica.fname_parse(845907)
        assert err.type == TypeError


def test_nii_convert():
    """
    """

    assert nii1 == meica.nii_convert(nii1)

    with pytest.raises(FileNotFoundError) as err:
        meica.nii_convert('empty.nii')
        assert err.type == FileNotFoundError

    assert 'nii_convert_test.nii.gz' == meica.nii_convert('nii_convert_test.nii')
    os.remove('nii_convert_test.nii.gz')

    try: # if 3dAFNItoNIFTI is available...
        assert 'sub-001_task-rest_run-02_echo-1_bold.nii.gz' == meica.nii_convert(afni1)
        os.remove('sub-001_task-rest_run-02_echo-1_bold.nii.gz')
    except: pass # unless AFNI is not in the current path


def test_format_inset():
    """
    """

    # NIFTI list dataset specification, string TE specification
    assert (nii2,
            'sub-001_task-rest_run-01_echo-123_bold') == meica.format_inset([nii1,nii2,nii3],
                                                                            ['12.2,24.6,30'])

    # NIFTI string dataset specification, list TE specification
    assert (nii2,
            'sub-001_task-rest_run-01_echo-123_bold') == meica.format_inset([','.join([nii1, nii2 ,nii3])],
                                                                            [12.2, 24.6, 30])

    # NIFTI shorthand dataset specification, string TE specification
    assert (nii2,
            'sub-001_task-rest_run-01_echo-123_bold') == meica.format_inset(['sub-001_task-rest_run-01_echo-[1,2,3]_bold.nii.gz'],
                                                                            ['12.2,24.6,30'])

    # AFNI shorthand dataset specification, list TE specification
    assert ('sub-001_task-rest_run-02_echo-2_bold+orig',
            'sub-001_task-rest_run-02_echo-123_bold') == meica.format_inset(['sub-001_task-rest_run-02_echo-[1,2,3]_bold+orig.BRIK.gz'],
                                                                            [12.2, 24.6, 30])

    # AFNI string dataset specification, list TE specification
    assert ('sub-001_task-rest_run-02_echo-2_bold+orig',
            'sub-001_task-rest_run-02_echo-123_bold') == meica.format_inset([','.join([afni1, afni2, afni3])],
                                                                            [12.2, 24.6, 30])

    # AFNI list dataset specification, string TE specification
    assert ('sub-001_task-rest_run-02_echo-2_bold+orig',
            'sub-001_task-rest_run-02_echo-123_bold') == meica.format_inset([afni1, afni2, afni3],
                                                                            ['12.2,24.6,30'])

    # catch differing number of datasets and echoes
    with pytest.raises(AssertionError) as err:
        meica.format_inset([','.join([nii1, nii2 ,nii3])],
                           [12.2, 24.6])
        assert err.type == AssertionError

    # catch incorrect echo input
    with pytest.raises(TypeError) as err:
        meica.format_tes(520885)
        assert err.type == TypeError

    # catch incorrect dataset input
    with pytest.raises(TypeError) as err:
        meica.format_dset(3258.27)
        assert err.type == TypeError


def test_find_CM():
    """
    """

    assert np.allclose([1.1718978881835938,
                        -41.01300048828125,
                        -46.293296813964844],
                        meica.find_CM('sub-001_T1w.nii.gz'))


def test_parser_opts():
    """
    """

    parser = meica.get_options(['-d', 'sub-001_task-rest_run-01_echo-[1,2,3]_bold.nii.gz',
                                '-e', '14.5,38.5,62.5'])

    # confirm required inputs are interpreted correctly from argparse
    assert parser.input_ds == ['sub-001_task-rest_run-01_echo-[1,2,3]_bold.nii.gz']
    assert parser.tes == ['14.5,38.5,62.5']

    # confirm that argparse complains about missing arguments
    with pytest.raises (SystemExit) as err:
        parser = meica.get_options(['-d', 'sub-001_task-rest_run-01_echo-[1,2,3]_bold.nii.gz'])
        assert err.type == SystemExit


# def test_gen_script():

    # os.environ['DYLD_FALLBACK_LIBRARY_PATH'] = '~/abin'
    # resdir = op.join(op.dirname(__file__),'resources')
    # fname = op.join(resdir,
    #                 '_meica_sub-001_task-rest_run-01_echo-123_bold.sh')
    # sel_opts = ['-d', 'sub-001_task-rest_run-01_echo-[1,2,3]_bold.nii.gz',
    #             '-e', '14.5,38.5,62.5',
    #             '-b', '4v',
    #             '-a', 'sub-001_T1w.nii.gz',
    #             '--fres=2', '--MNI', '--qwarp']

    # opts = meica.get_options(sel_opts)
    # with open(fname, 'r') as file:
    #     script = file.read()
    # script_list = meica.gen_script(opts)
    # "\n".join(script_list) == script
