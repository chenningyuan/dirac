

.. _molecule:

\**MOLECULE -- Specify the basis for the corresponding xyz file
===============================================================

This section only makes sense together with a xyz file.  The actual xyz file
only contains the coordinates.  Information about the atoms and the basis sets
is given in the inp file using this section.

The corresponding xyz file could look like this one::

  6
  Methanol, identifying the alcoholic proton as H1
  C      -0.000000010000       0.138569980000       0.355570700000   
  O      -0.000000010000       0.187935770000      -1.074466460000  
  H       0.882876920000      -0.383123830000       0.697839450000  
  H      -0.882876940000      -0.383123830000       0.697839450000  
  H      -0.000000010000       1.145042790000       0.750208830000  
  H1     -0.000000010000      -0.705300580000      -1.426986340000

The first line is the number of atoms.  The second line can be either empty or
contain a comment.  The following lines contain atom types and coordinates.
The first entry on each line is the name of the atom. Any string up to eight
characters will do, but if the name is present in the periodic table DIRAC will
recognize the charge and can find the corresponding basis set.

You can use the following keywords in this section:


\*BASIS
-------

Here the basis set information is specified.  .DEFAULT specifies the default
large component basis set for all atoms. This can be modified for specific atom
by the .SPECIAL keyword. For example, if we use the methanol.xyz file listed
above, the input::

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
section other basis sets can be specified in the traditional mol file.


\*CENTERS
---------

Here the physical properties of the centers in our molecule can be
specified. This is useful if the name of the center differs from the
IUPAC element name or in case we want to use nuclei with fractional
charges. This is controlled by the keyword nucleus that can be repeated
as often as is needed.

For the geometry file methanol.xyz we need to specify::

  *CENTERS
  .NUCLEUS
   H1 1.0


\*UNITS
-------

Here the units used for the coordinates in the xyz files can be entered.
If nothing is specified, angstrom are used (a real xyz file). If we
specify 'AU' we use atomic units instead::

  *UNITS
  AU

We can also specify our own unit, by typing the name and the conversion
factor when going from this new unit to atomic units::

  *UNITS
  nm 18.8972

Now all coordinates are entered in nanometer. If we add 'A' at the end,
we can enter the new unit in angstrom::

  *UNITS
  nm 10.0 A


\*SYMMETRY
----------

With this keyword we can control symmetry recognition. To turn symmetry
detection off and to force C1 symmetry use::

  *SYMMETRY
  .NOSYM
