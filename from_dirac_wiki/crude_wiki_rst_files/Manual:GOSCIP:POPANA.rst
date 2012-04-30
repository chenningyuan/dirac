Manual:GOSCIP:POPANA
====================

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

Contents
--------

-  `1 &POPANA -- Analyze the wave
   function <#.26POPANA_--_Analyze_the_wave_function>`_

   -  `1.1 THRESH <#THRESH>`_
   -  `1.2 DEGEN <#DEGEN>`_
   -  `1.3 SELPOP <#SELPOP>`_

&POPANA -- Analyze the wave function
====================================

`|image0| </dirac/index.php/Image:Dirac.png>`_

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

