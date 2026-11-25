Symmetry recognition in mol files
=================================

The test is meant to check the automatic detection of molecular symmetry in mol input files and calculates the nonrelativistic energy.

We assume that the basis sets are available in your standard path.
These sets can be obtained from the Dalton distribution and may later be included in the DIRAC distribution as well.

Note that a coarse DFT grid for LDA is used and the SCF calculation is stopped after 4 iterations.

Examples
--------
The following table contains few examples of
the symmetry detection from the xyz coordinates in mol files.


==============    =================
Molecule          Detected symmetry
==============    =================
CH4               T(d)  ->  D2
H2O               C(2v) ->  C2v 
SF6               O(h)  ->  D2h
trans CHF=CHF     C(2h) ->  C2h
CO_MgO_bridge     C(s)  ->  Cs
T                 T     ->  D2
Th                T(h)  ->  D2h
O                 O     ->  D2
I                 I     ->  D2
B12H12            I(h)  ->  D2h
==============    =================

``CO_MgO_bridge`` tests that centres are only symmetry-related if their
charges agree (to within a threshold of 1.0D-6), also for non-integer
charges, and that point charges (atom number 0) are distinguished by their
charge. It models CO on an oxide surface with embedding point charges:

* ``CO_MgO_bridge``: CO adsorbed upright at a bridge site of MgO(100), with the
  surface modelled by Mg(2+) and O(2-) point charges. The true symmetry is C(s)
  (the mirror plane containing the Mg-O bond). All point charges have atom
  number 0, and earlier versions found C(2v), imposing a mirror plane that
  swaps Mg(2+) and O(2-): the calculation ran without warning, but in the
  wrong field (energy off by 38 mEh).

The cubic and icosahedral cases (``T``, ``Th``, ``O``, ``I``, ``B12H12``) only run
the symmetry detection (``symdetect.inp``, no wave function).
``T``, ``Th``, ``O`` and ``I`` are a carbon atom surrounded by hydrogens on the orbit
of a general point, so no atom lies on a symmetry element. Earlier versions
reported T, O and I as T(d) and Th as O(h); in ``I`` the C5 axes were not found;
and icosahedral molecules were represented in C2 / C2h instead of D2 / D2h.

Few more examples can be added here to demonstrate (and to test) the wider variety of symmetry detections.
