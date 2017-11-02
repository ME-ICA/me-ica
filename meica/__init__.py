# emacs: -*- mode: python-mode; py-indent-offset: 4; tab-width: 4; indent-tabs-mode: nil -*-
# ex: set sts=4 ts=4 sw=4 et:
"""
MEICA: A Python package for multi-echo independent/principal component analysis.
"""

from duecredit import (due, Doi)

from .info import (
    __version__,
    __author__,
    __copyright__,
    __credits__,
    __license__,
    __maintainer__,
    __email__,
    __status__,
    __url__,
    __packagename__,
    __description__,
    __longdesc__
)

import warnings

# cmp is not used, so ignore nipype-generated warnings
warnings.filterwarnings('ignore', r'cmp not installed')

# Citation for the algorithm.
due.cite(Doi('10.1016/j.neuroimage.2011.12.028'),
         description='Introduces MEICA.',
         version=__version__, path='meica', cite_module=True)
due.cite(Doi('10.1073/pnas.1301725110'),
         description='Improves MEICA.',
         version=__version__, path='meica', cite_module=True)
