Manual:WAVE FUNCTION:RESOLVE
============================

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

Contents
--------

-  `1 \*RESOLVE -- Resolve open-shell
   states <#.2ARESOLVE_--_Resolve_open-shell_states>`_

   -  `1.1 .PRINT <#.PRINT>`_
   -  `1.2 .INTFLG <#.INTFLG>`_
   -  `1.3 .SCREEN <#.SCREEN>`_
   -  `1.4 .SCHEME <#.SCHEME>`_

\*RESOLVE -- Resolve open-shell states
======================================

The individual electronic states of an average-of-configurations
calculation are resolved by a full CI in the small open-shell
manifold(s) using the `GOSCIP module </dirac/index.php/Manual:GOSCIP>`_.

Note that even if this is a small CI diagonalization one needs in
principle to generate *all* AO integrals for the initial 4-index
transformation. The cost of calculation can be greatly reduced invoking
integral screening (controlled by `.SCREEN <#.SCREEN>`_) in the 4-index
transformation, possibly also by leaving out at least SS integrals
(using `.INTFLG <#.INTFLG>`_).

.PRINT
------

Print level.

*Default:*

::

    .PRINT
     0

.INTFLG
-------

Specify what two-electron integrals to include (see
`.INTFLG </dirac/index.php/Manual:HAMILTONIAN#.INTFLG>`_ under
`\*\*HAMILTONIAN </dirac/index.php/Manual:HAMILTONIAN>`_).

*Default:* `.INTFLG </dirac/index.php/Manual:HAMILTONIAN#.INTFLG>`_ from
`\*\*HAMILTONIAN </dirac/index.php/Manual:HAMILTONIAN>`_.

.SCREEN
-------

Screening threshold in the 4-index transformation.

A negative number deactivates screening.

*Default:*

::

    .SCREEN
     1.0D-14

.SCHEME
-------

Integral transformation algorithm.\ :sup:`[*please
`clarify </dirac/index.php/Clarify>`_*]`\ 

*Default:*

::

    .SCHEME
     4

