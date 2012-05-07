Manual:RELCCSD:CCENER
=====================

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

Contents
--------

-  `1 &CCENER -- Energy module <#.26CCENER_--_Energy_module>`_

   -  `1.1 DOMP2 <#DOMP2>`_
   -  `1.2 DOCCSD <#DOCCSD>`_
   -  `1.3 DOCCSDT <#DOCCSDT>`_
   -  `1.4 NTOL <#NTOL>`_
   -  `1.5 MAXIT <#MAXIT>`_
   -  `1.6 MAXDIM <#MAXDIM>`_
   -  `1.7 NOCCS <#NOCCS>`_
   -  `1.8 NOCCD <#NOCCD>`_

&CCENER -- Energy module
========================

DOMP2
-----

Calculate the MP2 energy.

*Default:*

::

    DOMP2=T

DOCCSD
------

Calculate the CCSD energy.

*Default:*

::

    DOCCSD=T

DOCCSDT
-------

Calculate the CCSD(T) energy. The program calculates three different
non-iterative triple corrections
\ :sup:``[1] <#_note-visscher:relccsd1>`_`\ .

*Default:*

::

    DOCCSDT=T

`|image0| </dirac/index.php/Image:Dirac.png>`_

**Advanced options**

NTOL
----

Convergence tolerance of the DIIS error vector (10\ :sup:`--NTOL`\ ).

*Default:*

::

    NTOL=12

MAXIT
-----

Maximum number of DIIS iterations to solve the CCSD equations.

*Default:*

::

    MAXIT=30

MAXDIM
------

Maximum dimension of the reduced space used in the DIIS extrapolation.

*Default:*

::

    MAXDIM=8

`|image1| </dirac/index.php/Image:HJAaJ.png>`_

**Programmers options**

NOCCS
-----

Ignore single excitations so that the CCD formalism is used. Only for
test purposes.

*Default:*

::

    NOCCS=F

NOCCD
-----

Ignore double excitations so that the CCS formalism is used. Only for
test purposes.

*Default:*

::

    NOCCD=F

