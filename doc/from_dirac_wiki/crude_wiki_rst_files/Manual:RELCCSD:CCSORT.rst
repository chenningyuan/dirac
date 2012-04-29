Manual:RELCCSD:CCSORT
=====================

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

Contents
--------

-  `1 &CCSORT -- Sorting module <#.26CCSORT_--_Sorting_module>`_

   -  `1.1 USEOE <#USEOE>`_
   -  `1.2 NORECMP <#NORECMP>`_

&CCSORT -- Sorting module
=========================

USEOE
-----

Use orbital energies from the HF calculation instead of the diagonal
elements of the reconstructed Fock matrix. Putting USEOE=F may solve
problems with low-symmetry open shell systems, in which a Kamers pair is
split over the occupied and virtual space, giving an orbital energy
difference of zero in the default definition.

*Default:*

::

    USEOE=T

NORECMP
-------

Do not recompute the Fock matrix for the determinant selected as
reference for the coupled cluster expansion. This is usually done as a
check on the numerical accuracy of the MO transformation. If you specify
both USEOE and NORECMP the reference Hamiltonian is given by the
original orbital energies (taken from the SCF program) and the
transformed 2-electron integrals. This makes it possible to apply the
molecular mean field approximation.

*Default:*

::

    NORECMP=F

