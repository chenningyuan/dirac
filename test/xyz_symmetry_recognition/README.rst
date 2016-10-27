Examples of xyz-based symmetry detection
========================================

We provide here (test/xyz_symmetry_recognition) collection of input tests aimed for the Dirac automatic symmetry recognition
of molecular structures given in the xyz-format.

Sources of xyz-inputs
---------------------

Molecules were obtained by a suitable software and 
their geometries (coordinates), reflecting the symmetry, was exported,
preferably in the xyz-format. 
If the GUI software does not save the geometry in xyz, try to save (export) in another
common format, and read it in an other software, capable to export the desired xyz file.

Very good choice for getting symmetrized xyz coordinates is the GUI of the ADF software.

Note that sometimes the resulting xyz-geometry file had to be manually controlled and coordinates 
corrected, if necessary, to keep the symmetry of the system.

List of symmetries
------------------

The hersym/FIND_PGROUP subroutine can handle these symmetries:
O(h); I(h); T(d); 
D(inf,h); D(nh); D(nd); D(n);
C(inf,v); C(nv); C(nh); C(n); C(s); C(i); C(1);
S(n)

Examples
--------

The following table contains examples of 
the symmetry detection in the xyz molecule inputs.

=========    ==========       ==========   ==============
xyz input    molecule         full group   represented as 
=========    ==========       ==========   ==============
molecule1    heterocyc        D(2h)        D2h
molecule2    Re2(CO)10        C(4v)        C2v
molecule3    Os(CO)5          C(2v)        C2v
molecule4    acetylene        D(oo,h)      D2h
molecule5    SF6              O(h)         D2h
molecule6    allene           D(2d)        C1 (wrong!) 
molecule7    NH3              C(3v)        Cs
molecule8    PF5              D(3h)        C2v
molecule9    Hs(CO)4          D(2d)        error exit 
=========    ==========       ==========   ===============

More examples should be added to this test to demostrate the full variety of symmetry detections.
