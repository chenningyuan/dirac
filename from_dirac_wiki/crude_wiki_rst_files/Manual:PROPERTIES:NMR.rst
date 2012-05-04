Manual:PROPERTIES:NMR
=====================

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

Contents
--------

-  `1 \*NMR -- Controls evaluation of NMR
   parameters <#.2ANMR_--_Controls_evaluation_of_NMR_parameters>`_

   -  `1.1 .LONDON <#.LONDON>`_
   -  `1.2 .GAUGEO <#.GAUGEO>`_
   -  `1.3 .GO\_ANG <#.GO_ANG>`_
   -  `1.4 .USECM <#.USECM>`_
   -  `1.5 .INTFLG <#.INTFLG>`_
   -  `1.6 .NOTWO <#.NOTWO>`_
   -  `1.7 .NOONEI <#.NOONEI>`_
   -  `1.8 .NOORTH <#.NOORTH>`_
   -  `1.9 .SYMCON <#.SYMCON>`_

\*NMR -- Controls evaluation of NMR parameters
==============================================

This section gives directives for the calculation of NMR parameters.

`|image0| </dirac/index.php/Image:Dirac.png>`_

**Advanced options**

.LONDON
-------

`|image1| </dirac/index.php/Image:Dirac.png>`_

**Only in the development version of DIRAC**

Activate calculations of magnetic properties (NMR shielding constants
and magnetizabilities with London atomic orbitals.

*Default:* Use conventional atomic orbitals.

.GAUGEO
-------

Reads in a user defined gauge origin (in bohr). Note that both the gauge
origin and dipole origin are changed in order to have consistency.

*Default for conventional atomic orbitals:*

::

    .GAUGEO
     0.0 0.0 0.0

.GO\_ANG
--------

Same as `.GAUGEO </dirac/index.php/Manual:PROPERTIES:NMR#.GAUGEO>`_ but
reads coordinates in angstrom.

.USECM
------

Use the center of mass as the gauge origin.

.INTFLG
-------

Specify what two-electron integrals to include in the two-electron
London contributions to the magnetic field property gradient (see
`.INTFLG </dirac/index.php/Manual:HAMILTONIAN#.INTFLG>`_ under
`\*\*HAMILTONIAN </dirac/index.php/Manual:HAMILTONIAN>`_).

*Default:* `.INTFLG </dirac/index.php/Manual:HAMILTONIAN#.INTFLG>`_ from
`\*\*HAMILTONIAN </dirac/index.php/Manual:HAMILTONIAN>`_.

.NOTWO
------

Do not calculate the two-electron London contributions for the magnetic
field property gradient when London atomic orbitals are used.

.NOONEI
-------

Do not calculate the H(0),T(B) reorthonormalization terms for the
magnetic field property gradient when London atomic orbitals are used.

.NOORTH
-------

Do not calculate the T(B),h(mK) reorthonormalization contributions for
the expectation value term when London atomic orbitals are used.

.SYMCON
-------

Employ the symmetric connection for reorthonormalization terms when
using London atomic orbitals.

*Default:* Use the natural connection.

