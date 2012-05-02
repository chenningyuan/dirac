Manual:ANALYZE:LOCALIZATION
===========================

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

Contents
--------

-  `1 \*LOCALIZATION -- Molecular orbital
   localization <#.2ALOCALIZATION_--_Molecular_orbital_localization>`_

   -  `1.1 .THRESH <#.THRESH>`_
   -  `1.2 .MAXITR <#.MAXITR>`_
   -  `1.3 .PRJLOC <#.PRJLOC>`_
   -  `1.4 .OWNBAS <#.OWNBAS>`_
   -  `1.5 .VECPRJ <#.VECPRJ>`_
   -  `1.6 .VECREF <#.VECREF>`_
   -  `1.7 .HESSUN <#.HESSUN>`_
   -  `1.8 .PRINT <#.PRINT>`_

\*LOCALIZATION -- Molecular orbital localization
================================================

Performs molecular orbitals localization using the Pipek-Mezey
criterion. The implementation uses an exponential parametrization and a
quasi-Newton minimization method instead of the Edminston-Ruedenberg
algorithm. For the moment, it only works with the *C*\ :sub:`1`\ 
symmetry and for the occupied orbitals.

.THRESH
-------

Convergence threshold.

*Default:*

::

    .THRESH
     1.0D-6

.MAXITR
-------

Maximum number of iterations.

*Default:*

::

    .MAXITR
     100

.PRJLOC
-------

Localization from projection analysis instead of Mulliken analysis.

::

    .PRJLOC
     1

*Default:*

::

    .PRJLOC
     0

.OWNBAS
-------

See `.OWNBAS </dirac/index.php/Manual:ANALYZE:PROJECTION#.OWNBAS>`_
under `\*PROJECTION </dirac/index.php/Manual:ANALYZE:PROJECTION>`_.

.VECPRJ
-------

See `.VECPRJ </dirac/index.php/Manual:ANALYZE:PROJECTION#.VECPRJ>`_
under `\*PROJECTION </dirac/index.php/Manual:ANALYZE:PROJECTION>`_.

.VECREF
-------

See `.VECREF </dirac/index.php/Manual:ANALYZE:PROJECTION#.VECREF>`_
under `\*PROJECTION </dirac/index.php/Manual:ANALYZE:PROJECTION>`_.

.HESSUN
-------

Hessian approximated by the unity matrix.

*Default:* Use the diagonal elements of the (diagonally dominant)
Hessian.

`|image0| </dirac/index.php/Image:HJAaJ.png>`_

**Programmers options**

.PRINT
------

Print level.

*Default:*

::

    .PRINT
     1

