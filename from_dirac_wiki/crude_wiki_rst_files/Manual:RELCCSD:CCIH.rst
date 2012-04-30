Manual:RELCCSD:CCIH
===================

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

Contents
--------

-  `1 &CCIH -- Intermediate Hamiltonian Fock space CC
   functionality <#.26CCIH_--_Intermediate_Hamiltonian_Fock_space_CC_functionality>`_

   -  `1.1 EHMIN <#EHMIN>`_
   -  `1.2 EHMAX <#EHMAX>`_
   -  `1.3 EPMIN <#EPMIN>`_
   -  `1.4 EPMAX <#EPMAX>`_

&CCIH -- Intermediate Hamiltonian Fock space CC functionality
=============================================================

EHMIN
-----

Lower boundary for the part of the *P*\ :sub:`*I*`\  space associated
with orbitals in which holes can be created (occupied orbitals for
sector 00).

*Default:*

::

    EHMIN = 1000.

EHMAX
-----

Upper boundary for the part of the *P*\ :sub:`*I*`\  space associated
with orbitals in which holes can be created (occupied orbitals for
sector 00).

*Default:*

::

    EHMAX = 100.

EPMIN
-----

Lower boundary for the part of the *P*\ :sub:`*I*`\  space associated
with orbitals in which particles can be created (virtual orbitals for
sector 00).

*Default:*

::

    EPMIN = -1.

EPMAX
-----

Upper boundary for the part of the *P*\ :sub:`*I*`\  space associated
with orbitals in which particles can be created (virtual orbitals for
sector 00).

*Default:*

::

    EPMAX = -1000.

