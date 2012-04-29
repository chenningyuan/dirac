Manual:Basis:XYZ
================

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

Contents
--------

-  `1 Introduction <#Introduction>`_
-  `2 Coordinates <#Coordinates>`_
-  `3 \*\*MOLECULE <#.2A.2AMOLECULE>`_

   -  `3.1 \*BASIS <#.2ABASIS>`_
   -  `3.2 \*CENTERS <#.2ACENTERS>`_
   -  `3.3 \*UNITS <#.2AUNITS>`_
   -  `3.4 \*SYMMETRY <#.2ASYMMETRY>`_

Introduction
------------

When a .xyz file is given to PAM it will automatically detect how to
handle this file. The actual .xyz file only contains the coordinates.
Information about the atoms and the basis sets is given in the .inp
file.

Coordinates
-----------

Coordinates are given in an .xyz file in the following way.

::

    3 # Number of atoms
    Water with 90 degree angles # Comment line 
    O 0.0   0.0   0.0
    H 1.0   1.0   0.0
    H 1.0 -1.0   0.0 

The first entry on each line is the name of the atom. Any string up to
eight characters will do, but if the name is present in the periodic
table Dirac will recognize the charge and can find the corresponding
basis set. In case you want to have a non-default basis set for a
specific atom you should define this in the \*\*MOLECULE block. E.g.
treat atom H1 different from the other Hs

::

    6
    Methanol, identifying the alcoholic proton as H1
    C      -0.000000010000       0.138569980000       0.355570700000   
    O      -0.000000010000       0.187935770000      -1.074466460000  
    H       0.882876920000      -0.383123830000       0.697839450000  
    H      -0.882876940000      -0.383123830000       0.697839450000  
    H      -0.000000010000       1.145042790000       0.750208830000  
    H1     -0.000000010000      -0.705300580000      -1.426986340000

\*\*MOLECULE
------------

The description of the atomic properties (basis, nuclear charge) is
given in the .INP file using the keyword \*\*MOLECULE when using
geometry (xyz) file input.

::

    **DIRAC
    Specification of the calculation.

::

    **MOLECULE
    Specification of the basis set and physical properties of nuclei (if different from the default).

In this section the following keywords are defined

\*BASIS
~~~~~~~

Here the basis set information is specified.

-  .DEFAULT

Specify the default `large
component </dirac/index.php/Manual:Basis:Large_component_basis_set>`_
basis set for all atoms. This can be modified for specific atom by the
.SPECIAL keyword. For example, if we use the methanol.xyz file listed
above, the input

::

    **MOLECULE
    *BASIS
    .DEFAULT
    cc-pVDZ
    .SPECIAL
    O BASIS cc-pVTZ
    .SPECIAL
    H1 BASIS cc-pVTZ
    *CENTERS
    .NUCLEUS
    H1 1.0
    *END OF INPUT

specifies a cc-pVDZ on carbon and the 3 hydrogens labeled H, while
oxygen and H1 will be treated with a cc-pVTZ basis. In this special
section other basis sets can be specified in the `same
way </dirac/index.php/Manual:Basis:Large_component_basis_set>`_ as in a
.MOL file

\*CENTERS
~~~~~~~~~

Here the physical properties of the centers in our molecule can be
specified. This is useful if the name of the center differs from the
IUPAC element name or in case we want to use nuclei with fractional
charges. This is controlled by the keyword nucleus that can be repeated
as often as is needed.

For the geometry file methanol.XYZ we need to specify

::

    *CENTERS
    .NUCLEUS
    H1 1.0

\*UNITS
~~~~~~~

Here the units used for the coordinates in the XYZ files can be entered.
If nothing is specified, angstroms are used (a real xyz file). If we
specify 'AU' we use atomic units instead. E.g.

::

    *UNITS
    AU

We can also specify our own unit, by typing the name and the conversion
factor when going from this new unit to atomic units.

::

    *UNITS
    nm 18.8972

Now all coordinates are entered in nanometer. If we add 'A' at the end,
we can enter the new unit in angstrom.

::

    *UNITS
    nm 10.0 A

\*SYMMETRY
~~~~~~~~~~

With this keyword you can control symmetry recognition. To turn symmetry
detection off and to force C1 symmetry use:

::

    *SYMMETRY
    .NOSYM

