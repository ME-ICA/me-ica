Usage
=====

ME-ICA minimally requires:

#. acquired echo times (in milliseconds), and
#. functional datasets equal to the number of acquired echoes.

But you can supply many other options, viewable with ``meica.py -h``.

Command line options
--------------------
.. argparse::
   :ref: meica.get_options(*)
   :prog: meica
   :nodefault:
   :nodefaultconst:

.. attention:: Make sure your datasets have slice timing information in the header.
   If not sure, specify a ``--tpattern`` option.
   Check AFNI documentation of `3dTshift`_ to see slice timing codes.

.. _3dTshift: http://afni.nimh.nih.gov/pub/dist/doc/program_help/3dTshift.html


.. tip:: FWHM smoothing is not recommended.
   tSNR boost is provided by optimal combination of echoes.
   For better overlap of 'blobs' across subjects, use non-linear standard space normalization with ``--qwarp``.
