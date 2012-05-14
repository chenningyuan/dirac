Manual:DIRRCI
=============

*Relativistic RASCI module* written by Lucas Visscher

This section allows for relativistic configuration interaction
calculations.

Note that when the namelist
`&GOSCIP </dirac/index.php/Manual:GOSCIP:GOSCIP>`_ is also present in
the input, `GOSCIP </dirac/index.php/Manual:GOSCIP>`_ will be called
instead of **DIRRCI**!

The input should be given in namelist form.


Manual:DIRRCI:CIROOT
====================

&CIROOT -- Select the state to converge on
==========================================

IREPNA
------

Abelian symmetry group of the desired state(s).

*Default: First abelian symmetry in the list*

NROOTS
------

Number of states to optimize on.

*Default:*

::

    NROOTS=1


**Advanced options**

ISTART
------

Start vector method.

COSCI start vectors:

::

    ISTART=1

Determinant with lowest eigenvalue:

::

    ISTART=2

First (reference) determinant:

::

    ISTART=3

If the COSCI start vectors are not available, then the default will be
2.

NSEL(NROOTS)
------------

Rank number of states to be optimized.

*Default:*

::

    NSEL=1,2,3,...

SELECT
------

Select wave functions on basis of largest overlap with the start wave
function.

*Default:*

::

    SELECT=F


Manual:DIRRCI:DIRECT
====================

&DIRECT -- Convergence control
==============================

CONVERE
-------

Convergence threshold for energy.

*Default:*

::

    CONVERR=1.0D-9

MAXITER
-------

Maximum number of direct CI iterations.

*Default:*

::

    MAXITER=10

**Advanced options**

CONVERR
-------

Convergence threshold for residual vector.

*Default:*

::

    CONVERR=1.0D-10

RESTART
-------

Restart on CI-vectors present in MRCFINV.

*Default:*

::

    RESTART=F

CPUMAX
------

Maximum amount of CPU-seconds to be used.

*Default:*

::

    CPUMAX=604800


Manual:DIRRCI:LEADDET
=====================

&LEADDET -- Analyze the CI wave function
========================================

**Advanced options**

GETDET
------

Get the list of dominant determinants.

*Default:*

::

    GETDET=T

COMIN
-----

Print contributions of determinants only if the square of the
coefficients is larger than COMIN.

*Default:*

::

    COMIN=0.1


Manual:DIRRCI:OPTIM
===================

&OPTIM -- Fine tuning of the algorithm
======================================

**Programmers options**

IGENEX
------

Write coupling coefficients to file (default):

::

    IGENEX=2

Calculate coupling coefficients when needed:

::

    IGENEX=1



Manual:DIRRCI:RASORB
====================

&RASORB -- Specify the type of CI and the active space
======================================================

NELEC
-----

Number of electrons (excluding frozen core electrons).

*Default:*

::

    NELEC=0

NRAS1
-----

Number of spinors in the RAS1 space for each abelian irrep.

*Default:*

::

    NRAS1=NSYMRP*0

NRAS2
-----

Number of spinors in the RAS2 space for each abelian irrep.

*Default:*

::

    NRAS2=NSYMRP*0

MAXH1
-----

Maximum number of holes in RAS1 spinors.

*Default:*

::

    MAXH1=0

MAXE3
-----

Maximum number of electrons in RAS3 spinors.

*Default:*

::

    MAXE3=0
