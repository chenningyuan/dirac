&GOSCIP -- Specify the CI space
===============================

NELEC
-----

Number of electrons (excluding frozen core electrons).

*Default:*

::

    NELEC=0

**NOTE: in the developer's version this line reads as**

NELACT
------

Number of electrons (excluding frozen core electrons).

*Default:*

::

    NELACT=0

**Programmers options**

IPRNT
-----

Print level.

*Default:*

::

    IPRNT=0


&POPANA -- Analyze the CI wave function
==========================================

**Advanced options**

THRESH
------

Print only determinants with coefficients higher than THRESH.

*Default:*

::

    THRESH=1.0D-3

DEGEN
-----

Threshold for when several states are considered to be degenerate.

*Default:*

::

    DEGEN=1.0D-10

SELPOP
------

Select only states with relative energies lower than SELPOP.

*Default:*

::

    SELPOP=1.0D2

