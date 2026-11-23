:orphan:


Restarting SCF calculations
===========================

A DIRAC calculation can be restarted at various stages depending on the type of
calculation. In general it should be enough to give the archive file of a
previous calculation as an argument to the ``pam`` script.  DIRAC will then
find and use the available data for restart.  Alternatively you can copy
restart files in using ``pam --put`` and if DIRAC finds them in the scratch
directory, it will automatically restart from them.  In general the program can
restart from

-  ``CHECKPOINT.h5`` Contains orbital coefficients in HDF5 format (machine-independent),
   and will be used for restarting the SCF procedure. Save with ``pam --outcmo``
   and restore with ``pam --incmo``.
-  ``DFPCMO`` Contains orbital coefficients in human-readable text format, and will be used for
   restarting the SCF procedure. Generated with the :ref:`GENERAL_.PCMOUT` keyword.
-  ``PAMFCK`` The Fock matrix. Will be used for restarting the SCF procedure.
-  ``AOPROPER`` Contains one-electron integrals.
-  ``X2CMAT`` Contains four- to two-component transformation matrices,
   and allows you to skip the transformation step of a two component
   calculation.
