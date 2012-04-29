Manual:DIRAC:OPTIMIZE
=====================

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

Contents
--------

-  `1 \*OPTIMIZE -- Geometry optimization
   directives <#.2AOPTIMIZE_--_Geometry_optimization_directives>`_

   -  `1.1 .NO SKIP <#.NO_SKIP>`_
   -  `1.2 .NUMGRA <#.NUMGRA>`_
   -  `1.3 .TWOGRD <#.TWOGRD>`_
   -  `1.4 .BAKER <#.BAKER>`_
   -  `1.5 .1STORD <#.1STORD>`_
   -  `1.6 Thresholds <#Thresholds>`_

      -  `1.6.1 .GRADIE <#.GRADIE>`_
      -  `1.6.2 .ENERGY <#.ENERGY>`_
      -  `1.6.3 .STEP T <#.STEP_T>`_

\*OPTIMIZE -- Geometry optimization directives
==============================================

This section controls the geometry optimization. The geometry
optimization algorithm is based on the one of
`DALTON <http://www.kjemi.uio.no/software/dalton/dalton.html>`_. The
directives are therefore the same except for these few changes:

-  No second-order algorithms are available since the molecular Hessian
   is not implemented.
-  A few new keywords have been introduced (see below).

 If `.PROPERTIES </dirac/index.php/Manual:DIRAC#.PROPERTIES>`_
(`.ANALYZE </dirac/index.php/Manual:DIRAC#.ANALYZE>`_) is specified in
the job input section together with
`.OPTIMIZE </dirac/index.php/Manual:DIRAC#.OPTIMIZE>`_, then the
property (analysis) module is called in each optimization iteration and
at the converged geometry.

In each iteration, the properties (analysis) requested under
`\*\*PROPERTIES </dirac/index.php/Manual:PROPERTIES>`_
(`\*\*ANALYZE </dirac/index.php/Manual:ANALYZE>`_) are calculated, and
at the converged geometry the properties (analysis) given in `\*\*PRP
F </dirac/index.php?title=Manual:PRP_F&action=edit>`_ (`\*\*ANA
F </dirac/index.php?title=Manual:ANA_F&action=edit>`_) are calculated.

.NO SKIP
--------

The `"trick" </dirac/index.php/Manual:PROPERTIES:MOLGRD#.TRICK>`_ is by
default activated in a geometry optimization, since when the current
geometry is far away from the equilibrium, e.g. the norm of the gradient
is, say, 1.0, then there is no need to calculate the LS and/or SS
two-electron gradient because they have a norm of, say 0.001.

This keywords forces the LS and SS two-electron gradient to always be
evaluated in all geometry iterations (depending on the integral flag).

.NUMGRA
-------

Force the use of numerical gradient in geometry optimization, e.g. in
Hartree-Fock calculations.

.TWOGRD
-------

*Default:* Include LL, SL, and SS integral contributions to the gradient
(1 = on; 0 = off).

::

    .TWOGRD 
     1 1 1

.BAKER
------

Baker's convergence criteria [J. Comp. Chem. 14(1993) 1085] will be
used.

.1STORD
-------

Default 1st order method will be used: BFGS update.

Thresholds
----------

Depending on convergence criterias for energies and gradients the values
of the thresholds may be automatically adjusted.

.GRADIE
~~~~~~~

Convergence threshold for the gradient. Default 1.0D-5.

.ENERGY
~~~~~~~

Convergence threshold for the energy. Default 1.0D-06.

.STEP T
~~~~~~~

Convergence threshold for the step. Default 0.10D-04.

