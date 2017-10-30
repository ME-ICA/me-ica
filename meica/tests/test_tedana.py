"""Tests for tedana."""

import os.path
import meica.tedana as tedana
import scripts.run_tedana as run_tedana

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
