"""Tests for tedana."""

import os.path
from meica.interfaces import tedana
from meica.cli import run_tedana

def test_basic_tedana():
    """
    A very simple test, to confirm that tedana creates output
    files.
    """

    parser = run_tedana.get_parser()
    options = parser.parse_args(['-d', '/home/neuro/data/zcat_ffd.nii.gz',
                                 '-e', '14.5', '38.5', '62.5'])
    tedana.main(options)
    assert os.path.isfile('comp_table.txt')
