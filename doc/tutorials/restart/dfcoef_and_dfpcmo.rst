:orphan:


What is the difference between CHECKPOINT.h5 and DFPCMO?
=========================================================

Both files hold the MO coefficients (and some other information) and can be
used to restart SCF calculations.

CHECKPOINT.h5 is always generated. To generate DFPCMO, you need the :ref:`GENERAL_.PCMOUT`
keyword.

To see the difference between the two files, run a calculation with::

  **GENERAL
  .PCMOUT

and save both files from the scratch directory::

  $ pam ... --outcmo --get "DFPCMO"

and then inspect both files. You will see that only DFPCMO is human-readable.
Both contain the **same** information.  The CHECKPOINT.h5 file is in HDF5 format
and is machine-independent, meaning that DIRAC can restart from a CHECKPOINT.h5
file generated on a different architecture. Since DIRAC 22, CHECKPOINT.h5 has
replaced the old binary DFCOEF file as the default checkpoint format. The
advantage of DFPCMO files is that you can read them yourself (try to understand
what the numbers mean) and use them in external scripts and programs.

You do not need any keywords for restarting: use ``pam --incmo`` to copy
CHECKPOINT.h5 to the scratch directory, or simply copy DFPCMO there,
and DIRAC will restart automatically.

Try it. It is important that you know how to restart. This will save you many
troubles and CPU hours.
