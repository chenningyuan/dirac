Manual:DIRRCI:CIROOT
====================

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

Contents
--------

-  `1 &CIROOT -- Select the state to converge
   on <#.26CIROOT_--_Select_the_state_to_converge_on>`_

   -  `1.1 IREPNA <#IREPNA>`_
   -  `1.2 NROOTS <#NROOTS>`_
   -  `1.3 ISTART <#ISTART>`_
   -  `1.4 NSEL(NROOTS) <#NSEL.28NROOTS.29>`_
   -  `1.5 SELECT <#SELECT>`_

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

`|image0| </dirac/index.php/Image:Dirac.png>`_

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

