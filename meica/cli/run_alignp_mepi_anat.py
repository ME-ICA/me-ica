# emacs: -*- mode: python-mode; py-indent-offset: 4; tab-width: 4; indent-tabs-mode: nil -*-
# ex: set sts=4 ts=4 sw=4 et:

import argparse
from meica import alignp_mepi_anat as ama


def get_parser():
    """
    Parses command line inputs

    Returns
    -------
    parser.parse_args() : argparse dic
    """
    parser = argparse.ArgumentParser()

    # Base processing options
    parser.add_argument('-e',
                        dest='tes',
                        nargs='+',
                        help='Echo times in ms. ex: -e 14.5 38.5 62.5',
                        required=True)
    parser.add_argument('-t',
                        dest='t2s',
                        help='T2* volume',
                        default='')
    parser.add_argument('-s',
                        dest='s0',
                        help='Skull-stripped S0 weighted volume, optional, for masking T2*',
                        default='')
    parser.add_argument('-a',
                        dest='anat',
                        help='Anatomical volume',
                        default='')
    parser.add_argument('-p',
                        dest='prefix',
                        help='Alignment matrix prefix',
                        default='')
    parser.add_argument('--cmass',
                        action='store_true',
                        dest='cmass',
                        help='Align cmass before main co-registration',
                        default=False)
    parser.add_argument('--autocmass',
                        action='store_false',
                        dest='autocmass',
                        help='Automatic cmass detection (default yes)',
                        default=True)
    parser.add_argument('--maxrot',
                        dest='maxrot',
                        help='Maximum rotation, default 30',
                        default='30')
    parser.add_argument('--maxshift',
                        dest='maxshf',
                        help='Maximum shift, default 30',
                        default='30')
    parser.add_argument('--maxscl',
                        dest='maxscl',
                        help='Maximum scale, default 1.01',
                        default='1.01')

    return parser


def main(argv=None):
    """Entry point"""
    options = get_parser().parse_args(argv)

    #Set up and cd into directory
    sl = []
    walignp_dirname = 'alignp.%s' % (options.prefix)
    sl.append('rm -rf %s ' % walignp_dirname)
    sl.append('mkdir %s ' % walignp_dirname)
    sl.append('cd %s ' % walignp_dirname)

    #Import datasets
    t2s_name, s0_name, anat_name = ama.import_datasets([options.t2s, options.s0, options.anat])

    #Filter and segment T2* volume to produce T2* base and weight volumes
    basevol, weightvol = ama.graywt(t2s_name, s0_name)

    #Run 3dAllineate
    allin_volume, _ = ama.allineate(anat_name, weightvol, basevol,
                                    options.prefix, options.maxrot,
                                    options.maxshf, options.maxscl,
                                    options.cmass)

    #For auto center-of-mass, check dice of COM and no-COM options
    ama.cmselect(basevol, allin_volume)

    #Run procedure
    ama.runproc(options.prefix, sl)


if __name__=='__main__':
    main()
