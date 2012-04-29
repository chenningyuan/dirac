Manual:RELCCSD:CCFOPR
=====================

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

Contents
--------

-  `1 &CCFOPR -- First-order
   properties <#.26CCFOPR_--_First-order_properties>`_

   -  `1.1 DOMP2G <#DOMP2G>`_
   -  `1.2 DOMP2GOLD <#DOMP2GOLD>`_
   -  `1.3 NEOPER <#NEOPER>`_
   -  `1.4 NAMEE <#NAMEE>`_
   -  `1.5 NTOL <#NTOL>`_
   -  `1.6 MAXIT <#MAXIT>`_
   -  `1.7 MAXDIM <#MAXDIM>`_

&CCFOPR -- First-order properties
=================================

Only MP2 first order properties are currently available. They are
computed by forming a relaxed density matrix, which is done (partly) in
AO basis to obtain all relaxation contributions. The properties that are
evaluated are those specified in the expectation value part of the
input.

DOMP2G
------

Calculate the MP2 gradient, that is, expectation values using the
relaxed one-electron density matrix. Note that this is not the molecular
gradient, but other first order properties as specified in the
expectation value part of the input.

*Default:*

::

    DOMP2G=F

`|image0| </dirac/index.php/Image:Dirac.png>`_

**Advanced options for the old formalism**

DOMP2GOLD
---------

Calculate the MP2 gradient via the old (approximate) MO formalism.

*Default:*

::

    DOMP2GOLD=F

NEOPER
------

Number of first order property operators |\\langle E \\rangle| (only
relevant for the old formalism).

*Default:*

::

    NEOPER=1

NAMEE
-----

Names of *E* operators in |image2|.

*Arguments:* Character (A8) (NAMEE(I),I = 1,NEOPER).

*Default:*

::

    NAMEE='ZDIPLEN'

NTOL
----

Convergence tolerance of the DIIS error vector (10\ :sup:`--NTOL`\ ).

*Default:*

::

    NTOL=12

MAXIT
-----

Maximum number of iterations to solve the Z-vector equations.

*Default:*

::

    MAXIT=30

MAXDIM
------

Maximum dimension of the reduced space.

*Default:*

::

    MAXDIM=8

