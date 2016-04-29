Symmetry detection
==================

Collection of input tests aimed for the automatic symmetry recognition
of molecular structures in the xyz-format.

Molecule were drawn by a suitable software. 
Then the geometry (coordinates), reflecting the symmetry, was exported,
preferably in the xyz-format. 
If the GUI software does not save the geometry in xyz, try to save (export) in another
common format, and read it in an other software, capable to export the desired xyz file.

Note that  the resulting xyz-geometry file had to be manually controlled and coordinates 
corrected, if necessary, to keep the symmetry of the system.

More examples can be added to demostrate the variety of symmetries
detected in the hersym.F/SYMADD subroutine.

xyz input    molecule    full group  represented as 
---------------------------------------------------
molecule1  aromatic napht.  D(2h)       D2h
molecule2  Re2(CO)10        C(4v)       C2v
molecule3  Os(CO)5          C(2v)       C2v
molecule4  acetylene        D(oo,h)     D2h
molecule5  SF6              O(h)        D2h
molecule6  allene           D(2d)       C1   <--- this is wrong !
molecule7  NH3              C(3v)       Cs
molecule8  PF5              D(3h)     

