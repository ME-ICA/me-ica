# emacs: -*- mode: python-mode; py-indent-offset: 4; tab-width: 4; indent-tabs-mode: nil -*-
# ex: set sts=4 ts=4 sw=4 et:
"""
MEICA: A Python package for multi-echo independent/principal component analysis.
"""

from .version import __version__
from .due import due, Doi, BibTeX

from .meica import afni_fname_parse
from .alignp_mepi_anat import dsprefix
from .t2smap import t2sadmap
from .tedana import do_svm
from .utils import cat2echos


__all__ = ['meica', 'tedana', 't2smap', 'alignp_mepi_anat', 'utils']

# Citation for the algorithm.
due.cite(Doi('10.1016/j.neuroimage.2011.12.028'),
         description='Introduces MEICA.',
         version=__version__, path='meica', cite_module=True)
due.cite(Doi('10.1073/pnas.1301725110'),
         description='Improves MEICA.',
         version=__version__, path='meica', cite_module=True)

# Cleanup
del due, Doi, BibTeX
del afni_fname_parse, dsprefix, t2sadmap, do_svm, cat2echos
