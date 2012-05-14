Restarting
==========

A DIRAC calculation can be restarted at various stages depending on the
type of calculation. In general it should be enough to give the archive
file of a previous calculation as an argument to the ``pam`` script.
DIRAC will then find and use the available data for restart. In general
the program can restart from

-  ``DFCOEF`` Contains orbital coefficients, and will be used for
   restarting the SCF procedure.
-  ``AOPROPER`` Contains one electron integrals.
-  ``BSSMAT`` Contains four- to two-component transformation matrices,
   and allows you to skip the make\_bss step of a two component
   calculation (if you DIRAC version supports this).

It is also possible to restart a correlated calculation, see for example
the manual section of `MOLTRA </dirac/index.php/Manual:MOLTRA>`_, and
also the tips on `how to run a coupled-cluster calculation in separate
steps </dirac/index.php/How_to_run_a_calculation_in_separate_steps>`_.

Restarting a Coupled-Cluster calculation
----------------------------------------

write me...
