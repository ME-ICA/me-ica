# -*- coding: utf-8 -*-
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""
Base module variables
"""

__version__ = '3.2.1-rc1-dev'
__author__ = ' meica developers'
__copyright__ = 'Copyright 2017, meica developers'
__credits__ = ['Elizabeth DuPre', 'Prantik Kundu', 'Ross Markello',
               'Taylor Salo', 'Kirstie Whitaker']
__license__ = 'LGPL 2.0'
__maintainer__ = 'Elizabeth DuPre'
__email__ = 'emd222@cornell.edu'
__status__ = 'Prototype'
__url__ = 'https://github.com/me-ica/me-ica'
__packagename__ = 'meica'
__description__ = ("Multi-echo Independent Components Analysis (ME-ICA) is a "
                   "functional magnetic resonance image pre-processing pipeline.")
__longdesc__ = ("To do.")

DOWNLOAD_URL = (
    'https://github.com/ME-ICA/{name}/archive/{ver}.tar.gz'.format(
        name=__packagename__, ver=__version__))

REQUIRES = [
    'numpy',
    'scikit-learn',
    'nilearn',
    'sklearn',
    'nibabel>=2.1.0',
    'pybids>=0.4.0',
    'nipype',
]

TESTS_REQUIRES = [
    "codecov",
    "pytest",
]

EXTRA_REQUIRES = {
    'doc': ['sphinx>=1.5.3', 'sphinx_rtd_theme', 'sphinx-argparse'],
    'tests': TESTS_REQUIRES,
    'duecredit': ['duecredit'],
}

# Enable a handle to install all extra dependencies at once
EXTRA_REQUIRES['all'] = [val for _, val in list(EXTRA_REQUIRES.items())]

# Package classifiers
CLASSIFIERS = [
    'Development Status :: 2 - Pre-Alpha',
    'Intended Audience :: Science/Research',
    'Topic :: Scientific/Engineering :: Image Recognition',
    'License :: OSI Approved :: LGPL License',
    'Programming Language :: Python :: 2.7',
    'Programming Language :: Python :: 3.6',
]
