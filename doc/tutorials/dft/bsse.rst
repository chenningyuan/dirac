
Basis set superposition error
=============================

The basis set superposition error (BSSE, also known as counterpoise correction) can be estimated by performing
subsystem calculations using ghost atoms which carry basis sets of the full system.

This tutorial describes how to use ghost atoms in DIRAC and how to avoid
certain pitfalls when working with ghost centers.

Treating BSSE in xyz-files
~~~~~~~~~~~~~~~~~~~~~~~~~~
The xyz-reader offers a simple way for dealing with counterpoise procedure.
Let us have a look at this example. The Be2.xyz molecule input file labels the second Berylium atom 
as 'ghost' (what means zero charge) ::

 Beryllium dimer
 Be   0.0 0.0 0.0
 BeGh 0.0 0.0 4.0

The corresponding Dirac input file has to contain the reference to the ghost nucleus with two numbers. 
The first specifies the atom's nuclear charge that is used in the calculation (here it is zero), 
and the second one is the proton number used for searching in basis set files::

 **DIRAC
 .TITLE
  Be dimer counterpoise correction
 .WAVE FUNCTION
 **WAVE FUNCTION
 .SCF
 *SCF
 .CLOSED SHELL
  4
 .EVCCNV
 2.0d-07 2.0d-07
 **MOLECULE
 *CENTERS
 # 2 numbers - nuclear charge for ghost atom and the proton number
 .NUCLEUS
 BeGh 0.0 4.0
 *BASIS
 .DEFAULT 
 cc-pVDZ
 **END OF

For more about specifying xyz-files see the :ref:`molecule_using_xyz` part of the manual.

Using basis set file in mol-files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
If you are using the mol-file input with reading of basis set(s) from the file(s),
you can use this simple scheme::

 INTGRL
  Ne2   neon dimer
  Using aug-cc-pVQZ-basis
 C   2    2  X  Y
      0.0     1
 Ne  0.0 0.0  -2.55
 LARGE BASIS Q=10.0  aug-cc-pVQZ
     10.0      1
 Ne2  0.0 0.0   2.55
 LARGE BASIS aug-cc-pVQZ
 FINISH

The first neon atom is considered as the 'ghost' atom with the charge of 0.0 in the mol-file, but also with the "Q=10.0" 
parameter before the basis set file name (aug-cc-pVQZ). The "Q=..." value represents the real charge of that atom so
that the program can find it in user's specified basis set file.

Explicit giving of ghost atom basis set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Last (and least comfortable) way is to type the ghost atom basis set explicitly. This is when there is no
basis set file at disposal, and the user has to write basis set exponents (and contraction coefficients)
directly into the mol-file.

As to the replacement of a real atom by its ghost atom, it is sufficient to
replace the atom's charge by zero charge. It is recommended to also change the
name of the atom. 

For instance, this input for a ghost atom with a contracted
oxygen aug-cc-pVDZ basis could be given by the following molecule (.mol) file::

          0.    1
  Gh      .0000000000        0.0000000000        -.2249058930
  LARGE EXPLICIT  3    1    1    1
  # s functions
  f  10    4
    11720.0000000  0.00071000 -0.00016000  0.00000000  0.00000000
     1759.0000000  0.00547000 -0.00126300  0.00000000  0.00000000
      400.8000000  0.02783700 -0.00626700  0.00000000  0.00000000
      113.7000000  0.10480000 -0.02571600  0.00000000  0.00000000
       37.0300000  0.28306200 -0.07092400  0.00000000  0.00000000
       13.2700000  0.44871900 -0.16541100  0.00000000  0.00000000
        5.0250000  0.27095200 -0.11695500  0.00000000  0.00000000
        1.0130000  0.01545800  0.55736800  0.00000000  0.00000000
        0.3023000 -0.00258500  0.57275900  1.00000000  0.00000000
        0.0789600  0.00000000  0.00000000  0.00000000  1.00000000
  # p functions
  f   5    3
       17.7000000  0.04301800  0.00000000  0.00000000
        3.8540000  0.22891300  0.00000000  0.00000000
        1.0460000  0.50872800  0.00000000  0.00000000
        0.2753000  0.46053100  1.00000000  0.00000000
        0.0685600  0.00000000  0.00000000  1.00000000
  # d functions
  f   2    2
        1.1850000  1.00000000  0.00000000
        0.3320000  0.00000000  1.00000000

For detailed know-how about mol files see :ref:`molecule_using_mol`.

DFT grid
--------

When you replace an atom by a ghost atom the DFT grid will change. In order to
use the same grid for both the full system and the subsystem with ghost atoms,
you can export the DFT grid of the full system by copying the file
"numerical\_grid" from the scratch directory, and import this grid back into
the calculation with ghost atoms::

  **GRID
  .IMPORT
  numerical_grid

Make sure to copy the grid file from/to the working directory by using::

  pam ... --get=numerical_grid ...

and::

  pam ... --copy=numerical_grid ...


Symmetry recognition
--------------------

If you leave the symmetry text in the molecule specification blank, DIRAC will
try to evaluate the symmetry of the system automatically.  Note that the
symmetry recognition uses the atom types and coordinates but not the respective
basis sets. This means that if you use ghost atoms with different basis sets,
you should give the symmetry explicitly or check the symmetry recognized by
DIRAC.
