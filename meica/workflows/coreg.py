from nipype.interfaces import fsl
import nipype.interfaces.afni as afni
from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu
from nipype.interfaces.utility import Function


def t2s_coreg_wf(name='t2s_coreg_wf'):
    """
    This workflow .
    Outputs are .
    .. workflow::
        :graph2use: colored
        :simple_form: yes
        from meica.workflows.coreg import t2s_coreg_wf

    Parameters
    ----------
    name: str
        Workflow name (default: t2s_coreg_wf)

    Inputs
    ----------
    t1w
        List of T1w structural images
    t2s
        List of T2* functional images

    Outputs
    ----------
    coreg_params
        Affine transform

    """

    workflow = pe.Workflow(name='t2s_coreg_wf')

    inputnode = pe.Node(niu.IdentityInterface(fields=['t2svol', 'anat']),
                        name='inputnode')

    outputnode = pe.Node(niu.IdentityInterface(fields=['coreg_params']),
                         name='outputnode')

    get_thr = pe.Node(fsl.ImageStats(op_string='-P 50'),
                      name='get_thr')

    def format_expr(val):
        """
        Generates string for use as `expr`
        input in afni.Calc()

        Parameters
        ----------
        val: float
            Threshold generated from fsl.ImageStats()

        Outputs
        ----------
        expr_string
            Expression to be applyed with afni.Calc()
        """
        expr_string = 'a*isnegative(a-2*{})'.format(val)
        return expr_string

    fmt_expr = pe.Node(name='fmt_expr',
                       interface=Function(input_names=['val'],
                                          output_names=['expr_string'],
                                          function=format_expr))

    apply_thr = pe.Node(afni.Calc(), name='apply_thr')

    t1_seg = pe.Node(fsl.FAST(use_priors=True,
                              probability_maps=True), name='t1_seg')

    align = pe.Node(afni.Allineate(out_file='mepi_al.nii.gz',
                                   out_matrix='mepi_al_mat.1D',
                                   source_automask=2,
                                   warp_type='affine_general',
                                   args='-weight_frac 1.0 -lpc'
                                        '-maxshf 30 -maxrot 30'
                                        '-maxscl 1.01'),
                    name='align')

    workflow.connect([
                    (inputnode, get_thr, [('t2svol', 'in_file')]),
                    (inputnode, align, [('anat', 'in_file'),
                                        ('anat', 'master')]),
                    (get_thr, fmt_expr, [('out_stat', 'val')]),
                    (inputnode, apply_thr, [('t2svol', 'in_file_a')]),
                    (fmt_expr, apply_thr, [('expr_string', 'expr')]),
                    (apply_thr, t1_seg, [('out_file', 'in_files')]),
                    (apply_thr, align, [('out_file', 'reference')]),
                    (t1_seg, align, [('tissue_class_map', 'weight_file')]),
                    (align, outputnode, [('matrix', 'coreg_params')])
                    ])

    workflow.write_graph(graph2use='colored', simple_form=True)

    return workflow


t2s_coreg_wf()
