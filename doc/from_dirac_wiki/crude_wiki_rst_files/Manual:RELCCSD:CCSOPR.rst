Manual:RELCCSD:CCSOPR
=====================

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

Contents
--------

-  `1 &CCSOPR -- Second-order
   properties <#.26CCSOPR_--_Second-order_properties>`_

   -  `1.1 NAOPER <#NAOPER>`_
   -  `1.2 NBOPER <#NBOPER>`_
   -  `1.3 NFREQ <#NFREQ>`_
   -  `1.4 NAMEA <#NAMEA>`_
   -  `1.5 NAMEB <#NAMEB>`_
   -  `1.6 EFREQ <#EFREQ>`_
   -  `1.7 NTOL <#NTOL>`_
   -  `1.8 MAXIT <#MAXIT>`_
   -  `1.9 MAXDIM <#MAXDIM>`_

&CCSOPR -- Second-order properties
==================================

Only the SCF wave function can be used at the moment (RPA). The names of
the property integrals should first be defined in the
`.PRPTRA </dirac/index.php/Manual:MOLTRA#.PRPTRA>`_ section of
`\*\*MOLTRA </dirac/index.php/Manual:MOLTRA>`_.

NAOPER
------

Number of A operators in |\\langle\\langle A ; B
\\rangle\\rangle\_\\omega|.

*Default:*

::

    NAOPER=1

NBOPER
------

Number of B operators in |image1|.

*Default:*

::

    NBOPER=1

NFREQ
-----

Number of frequencies in |image2|.

*Default:*

::

    NFREQ=1

NAMEA
-----

Names of A operators in |image3|.

*Arguments:* Character (A8) (NAMEA(I),I = 1,NAOPER).

*Default:*

::

    NAMEA='ZDIPLEN'

NAMEB
-----

Names of B operators in |image4|.

*Arguments:* Character (A8) (NAMEB(I),I = 1,NBOPER).

*Default:*

::

    NAMEB='ZDIPLEN'

EFREQ
-----

Values of frequencies in |image5|.

*Arguments:* Real (EFREQ(I),I = 1,NFREQ).

*Default:*

::

    EFREQ=0.0

`|image6| </dirac/index.php/Image:Dirac.png>`_

**Advanced options**

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

