Examples of xyz-based symmetry detection
========================================

We provide here (see test/xyz_symmetry_recognition) collection of input tests aimed for the Dirac automatic symmetry recognition
of molecular structures given in the xyz-coordinates format.

Sources of xyz-inputs
---------------------

Molecules were obtained by a suitable software and 
their geometries (coordinates), reflecting the symmetry, was exported,
preferably in the xyz-format. 
If the GUI software does not save the geometry in xyz, try to save (export) in another
common format, and read it in an other software, capable to export the desired xyz file.
Note that sometimes the resulting xyz-geometry file had to be manually controlled and coordinates 
corrected, if necessary, to keep the symmetry of the system.

Very good choice for getting symmetry imposed xyz coordinates is the GUI of the commercial ADF software.

List of symmetries
------------------

The DIRAC's *hersym.F/FIND_PGROUP* subroutine can handle the following point (full) groups:

* O(h) 
* I(h) 
* T(d)
* D(inf,h) 
* D(nh)
* D(nd)
* D(n)
* C(inf,v)
* C(nv)
* C(nh)
* C(n)
* C(s)
* C(i)
* C(1)
* S(n)

The detected full group is then turned into lower computational group,
which represents the found point group.


Examples
--------

The following table contains examples of 
the symmetry detection from the xyz coordinates.

==========   =====================     ==========    ==============
xyz input    molecule                  full group    represented as 
==========   =====================     ==========    ==============
molecule01   heterocycycle             D(2h)         D2h
molecule02   Re2(CO)10                 C(4v)         C2v
molecule03   Os(CO)5                   C(2v)         C2v
molecule04   acetylene                 D(oo,h)       D2h
molecule05   SF6                       O(h)          D2h
molecule06   allene, S(4)              D(2d)         D2
molecule07   NH3                       C(3v)         Cs
molecule08   PF5                       D(3h)         C2v
molecule09   Hs(CO)4, D(2d)            D(2) (!)      error exit! 
molecule10   ethane                    C(2h)         C2h
molecule11   dodecahedrane, I(h)       C(2h) (!)     C2h
molecule12   dodecahedrane, D(5d)      D(5d)         C2h
molecule13   benzene, D(6h)            D(6h)         D2h 
==========   =====================     ==========    ==============

More examples have to be be added to this test to demostrate (and test) the complete variety of symmetry detections.
