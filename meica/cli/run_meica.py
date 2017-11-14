# emacs: -*- mode: python-mode; py-indent-offset: 4; tab-width: 4; indent-tabs-mode: nil -*-
# ex: set sts=4 ts=4 sw=4 et:

import argparse
import meica.workflows.meica as meica


def get_parser():
    """
    Parses command line inputs

    Returns
    -------
    parser.parse_args() : argparse dic
    """

    # Configure options and help dialog
    parser = argparse.ArgumentParser()

    # Base processing options
    parser.add_argument('-e',
                        dest='tes',
                        nargs='+',
                        help='Echo times in ms. ex: -e 14.5 38.5 62.5',
                        required=True)
    parser.add_argument('-d',
                        dest='input_ds',
                        nargs='+',
                        help='Input datasets. ex: -d RESTe[123].nii.gz',
                        required=True)
    parser.add_argument('-a',
                        dest='anat',
                        help='(Optional) anatomical dataset. ' +
                             'ex: -a mprage.nii.gz',
                        default='')
    parser.add_argument('-b',
                        dest='basetime',
                        help='Time to steady-state equilibration in ' +
                             'seconds(s) or volumes(v). Default 0. ex: -b 4v',
                        default='0')
    parser.add_argument('--MNI',
                        dest='mni',
                        action='store_true',
                        help='Warp to MNI standard space.',
                        default=False)
    parser.add_argument('--strict',
                        dest='strict',
                        action='store_true',
                        help='Use strict component selection, suitable with' +
                             ' large-voxel, high-SNR data',
                        default=False)

    # Extended options for processing
    extopts = parser.add_argument_group('Additional processing options')
    extopts.add_argument('--qwarp',
                         dest='qwarp',
                         action='store_true',
                         help='Nonlinear warp to standard space using QWarp,' +
                              ' requires --MNI or --space).',
                         default=False)
    extopts.add_argument('--native',
                         dest='native',
                         action='store_true',
                         help='Output native space results in addition to ' +
                              'standard space results.',
                         default=False)
    extopts.add_argument('--space',
                         dest='space',
                         help='Path to specific standard space template for ' +
                              'affine anatomical normalization.',
                         default=False)
    extopts.add_argument('--fres',
                         dest='fres',
                         help='Specify functional voxel dimensions in mm ' +
                              '(iso.) for resampling during preprocessing.' +
                              'ex: --fres=2.5',
                         default=False)
    extopts.add_argument('--no_skullstrip',
                         action='store_true',
                         dest='no_skullstrip',
                         help='Anat is intensity-normalized and ' +
                              'skull-stripped (for use with -a flag).',
                         default=False)
    extopts.add_argument('--no_despike',
                         action='store_true',
                         dest='no_despike',
                         help='Do not de-spike functional data. ' +
                              'Default is to de-spike, recommended.',
                         default=False)
    extopts.add_argument('--no_axialize',
                         action='store_true',
                         dest='no_axialize',
                         help='Do not re-write dataset in axial-first order.' +
                              ' Default is to axialize, recommended.',
                         default=False)
    extopts.add_argument('--mask_mode',
                         dest='mask_mode',
                         help='Mask functional with help from anatomical or' +
                              ' standard space images.' +
                              ' Options: "anat" or "template".',
                         default='func')
    extopts.add_argument('--coreg_mode',
                         dest='coreg_mode',
                         help='Coregistration with Local Pearson and T2* weights '+
                              '(default), or use align_epi_anat.py (edge method).'+
                              'Options: "lp-t2s" or "aea"',
                         default='lp-t2s')
    extopts.add_argument('--smooth',
                         dest='FWHM',
                         help='FWHM smoothing with 3dBlurInMask. Default none. ' +
                              'ex: --smooth 3mm ',
                         default='0mm')
    extopts.add_argument('--align_base',
                         dest='align_base',
                         help='Explicitly specify base dataset for volume ' +
                              'registration',
                         default='')
    extopts.add_argument('--TR',
                         dest='TR',
                         help='TR. Read by default from dataset header',
                         default='')
    extopts.add_argument('--tpattern',
                         dest='tpattern',
                         help='Slice timing (i.e. alt+z, see 3dTshift -help).' +
                              ' Default from header.',
                         default='')
    extopts.add_argument('--align_args',
                         dest='align_args',
                         help='Additional arguments to anatomical-functional' +
                              ' co-registration routine',
                         default='')
    extopts.add_argument('--ted_args',
                         dest='ted_args',
                         help='Additional arguments to ' +
                              'TE-dependence analysis',
                         default='')

    # Additional, extended preprocessing options
    #  no help provided, caveat emptor
    extopts.add_argument('--select_only',
                         dest='select_only',
                         action='store_true',
                         help=argparse.SUPPRESS,
                         default=False)
    extopts.add_argument('--tedica_only',
                         dest='tedica_only',
                         action='store_true',
                         help=argparse.SUPPRESS,
                         default=False)
    extopts.add_argument('--export_only',
                         dest='export_only',
                         action='store_true',
                         help=argparse.SUPPRESS,
                         default=False)
    extopts.add_argument('--daw',
                         dest='daw',
                         help=argparse.SUPPRESS,
                         default='10')
    extopts.add_argument('--tlrc',
                         dest='space',
                         help=argparse.SUPPRESS,
                         default=False)  # For backwards compat.
    extopts.add_argument('--highpass',
                         dest='highpass',
                         help=argparse.SUPPRESS,
                         default=0.0)
    extopts.add_argument('--detrend',
                         dest='detrend',
                         help=argparse.SUPPRESS,
                         default=0.)
    extopts.add_argument('--initcost',
                         dest='initcost',
                         help=argparse.SUPPRESS,
                         default='tanh')
    extopts.add_argument('--finalcost',
                         dest='finalcost',
                         help=argparse.SUPPRESS,
                         default='tanh')
    extopts.add_argument('--sourceTEs',
                         dest='sourceTEs',
                         help=argparse.SUPPRESS,
                         default='-1')
    parser.add_argument_group(extopts)

    # Extended options for running
    runopts = parser.add_argument_group('Run options')
    runopts.add_argument('--prefix',
                         dest='prefix',
                         help='Prefix for final ME-ICA output datasets.',
                         default='')
    runopts.add_argument('--cpus',
                         dest='cpus',
                         help='Maximum number of CPUs (OpenMP threads) to use. ' +
                         'Default 2.',
                         default='2')
    runopts.add_argument('--label',
                         dest='label',
                         help='Label to tag ME-ICA analysis folder.',
                         default='')
    runopts.add_argument('--test_proc',
                         action='store_true',
                         dest='test_proc',
                         help='Align, preprocess 1 dataset then exit.',
                         default=False)
    runopts.add_argument('--pp_only',
                         action='store_true',
                         dest='preproc_only',
                         help='Preprocess only, then exit.',
                         default=False)
    runopts.add_argument('--keep_int',
                         action='store_true',
                         dest='keep_int',
                         help='Keep preprocessing intermediates. ' +
                              'Default delete.',
                         default=False)
    runopts.add_argument('--RESUME',
                         dest='resume',
                         action='store_true',
                         help='Attempt to resume from normalization step ' +
                              'onwards, overwriting existing')
    runopts.add_argument('--OVERWRITE',
                         dest='overwrite',
                         action='store_true',
                         help='Overwrite existing meica directory.',
                         default=False)
    parser.add_argument_group(runopts)

    return parser


def main(argv=None):
    """Entry point"""
    options = get_parser().parse_args(argv)

    if options.debug:
        logger.setLevel(logging.DEBUG)

    script_list = meica.gen_script(options)

    setname = meica.get_setname(options.input_ds, options.tes)
    if options.label != '':
        setname = setname + options.label

    with open('_meica_{}.sh'.format(setname), 'w') as file_object:
        file_object.write('\n'.join(script_list))


if __name__ == '__main__':
    main()
