Manual:PROPERTIES:EXPECTATION VALUE
===================================

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

Contents
--------

-  `1 \*EXPECTATION VALUE -- Evaluate expectation
   values <#.2AEXPECTATION_VALUE_--_Evaluate_expectation_values>`_

   -  `1.1 .OPERATOR <#.OPERATOR>`_
   -  `1.2 .ORBANA <#.ORBANA>`_
   -  `1.3 .PROJEC <#.PROJEC>`_
   -  `1.4 .PRINT <#.PRINT>`_

\*EXPECTATION VALUE -- Evaluate expectation values
==================================================

This section gives directives for the calculation of expectation values.

.OPERATOR
---------

Specification of general one-electron operators.

See the `specification of one-electron
operators </dirac/index.php/Manual:One-electron_operators>`_ section for
more information and explicit examples.

`|image0| </dirac/index.php/Image:Dirac.png>`_

**Advanced options**

.ORBANA
-------

Analyze contributions to an expectation value from the individual
orbitals.

.PROJEC
-------

Decompose an expectation value in contribution from fragment orbitals.
This action presupposes that a projection analysis has been performed,
see the `projection
analysis </dirac/index.php/Manual:ANALYZE:PROJECTION>`_ section for
details.

`|image1| </dirac/index.php/Image:HJAaJ.png>`_

**Programmers options**

.PRINT
------

Print level.

*Default:*

::

    .PRINT
     0

