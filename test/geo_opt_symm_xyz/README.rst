Geometry optimization tests
===========================

Collection of geometry optimizations tests of symmetric molecules 
with xyz molecule inputs.

Testsa are aimed at checking the symmetry detection and proper geometry displacements
of simple symmetric molecules.

Resulting converged geometries must be of the given symmetry.

CH4
---

Has the T(d) point group symmetry.
Sets Cs computational symmetry.

NH3
---

Has the C(3v) point group symmetry.
Sets Cs computational symmetry.

Allene
------

Has the D(2d) point group symmetry.
Molecule can not detect symmetry and during optimization 
it gets wrong coordinates resulting in SCF nonconvergence.


SF5Cl
-----

Has the C(4v) point group symmetry.
The symmetry detection (C2v) ends with error - Nuclei too close.


SF6
-----

Has the O(h) point group symmetry.
Symmetry detection ends with error - Nuclei too close.


