Manual:PROPERTIES:MOLGRD
========================

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

Contents
--------

-  `1 \*MOLGRD -- Controls evaluation of the molecular
   gradient <#.2AMOLGRD_--_Controls_evaluation_of_the_molecular_gradient>`_

   -  `1.1 .TRICK <#.TRICK>`_
   -  `1.2 .INTFLG <#.INTFLG>`_
   -  `1.3 .SCREEN <#.SCREEN>`_
   -  `1.4 .PRINT <#.PRINT>`_
   -  `1.5 .NUMGRA <#.NUMGRA>`_

\*MOLGRD -- Controls evaluation of the molecular gradient
=========================================================

The section gives directives for control of a single point evaluation of
the molecular gradient activated with
`.MOLGRD </dirac/index.php/Manual:PROPERTIES#.MOLGRD>`_ under
`\*\*PROPERTIES </dirac/index.php/Manual:PROPERTIES>`_.

The molecular gradient module is usually called from within the geometry
optimization module.

`|image0| </dirac/index.php/Image:Dirac.png>`_

**Advanced options**

.TRICK
------

If activated the calculation of the SL and SS two-electron integral
contributions to the gradient is skipped if their contribution is
estimated to be small. An "empirical" estimate of the norms of the LS
and SS two-electron gradients based on the norm of the LL two-electron
gradient. If deactivated all contributions to gradient are calculated
and if activated then only "necessary" contributions are calculated.

*Default:* Deactivated.

.INTFLG
-------

Specify what two-electron integrals to include (see
`.INTFLG </dirac/index.php/Manual:HAMILTONIAN#.INTFLG>`_ under
`\*\*HAMILTONIAN </dirac/index.php/Manual:HAMILTONIAN>`_).

*Default:* `.INTFLG </dirac/index.php/Manual:HAMILTONIAN#.INTFLG>`_ from
`\*\*HAMILTONIAN </dirac/index.php/Manual:HAMILTONIAN>`_.

.SCREEN
-------

Screening threshold.

*Default:*

::

    .SCREEN
     1.0D-17

Which means that screening is turned off by default. This is due to an
enormous demand of memory.

`|image1| </dirac/index.php/Image:HJAaJ.png>`_

**Programmers options**

.PRINT
------

Print level.

A print level of less than 3 gives only the total gradient. Print levels
from 3+ gives individual contributions to the gradient (kinetic energy
gradient, nuclear attraction gradient etc.). Print levels of 5+ print
some matrices, and 10+ give massive output!

*Default:* Print level from
`\*\*GENERAL </dirac/index.php/Manual:GENERAL>`_.

.NUMGRA
-------

Calculate the molecular gradient using numerical differentiation.

*Default:* Analytical evaluation if implemented.

