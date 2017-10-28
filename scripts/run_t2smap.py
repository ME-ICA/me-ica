# emacs: -*- mode: python-mode; py-indent-offset: 4; tab-width: 4; indent-tabs-mode: nil -*-
# ex: set sts=4 ts=4 sw=4 et:
import argparse

from meica import t2smap


def get_parser():
    """
    Parses command line inputs for t2smap
    Returns
    -------
    parser.parse_args() : argparse dic
    """
    parser = argparse.ArgumentParser()
    parser.add_option('-d',
                      nargs='+',
                      dest='data',
                      help='Spatially Concatenated Multi-Echo Dataset',
                      required=True)
    parser.add_option('-e',
                      nargs='+',
                      dest='tes',
                      help='Echo times (in ms) ex: 15,39,63',
                      required=True)
    parser.add_option('-c',
                      dest='combmode',
                      help='Combination scheme for TEs: t2s (Posse 1999),ste(Poser,2006 default)',
                      default='ste')
    parser.add_option('-l',
                      dest='label',
                      help='Optional label to tag output files with',
                      default=None)
    return parser


def main(argv=None):
    parser = get_parser()
    options = parser.parse_args(argv)
    t2smap.main(options)


if __name__=='__main__':
    main()
