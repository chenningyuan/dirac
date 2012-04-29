Manual:Basis:Atomic coordinates
===============================

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

Contents
--------

-  `1 Atomic coordinates <#Atomic_coordinates>`_

   -  `1.1 Proton number <#Proton_number>`_
   -  `1.2 Number of unique atoms <#Number_of_unique_atoms>`_
   -  `1.3 Nuclear exponent <#Nuclear_exponent>`_
   -  `1.4 Atom name <#Atom_name>`_
   -  `1.5 Cartesian coordinates <#Cartesian_coordinates>`_

Atomic coordinates
==================

Note that only the top six lines are shown in the input examples below.

Proton number
-------------

::

    DIRAC
    Water
    aug-cc-pVDZ basis (note that a tight p has NOT been added)
    C   2    2  X  Y
            8.    1
    O      .0000000000        0.0000000000        -.2249058930
    ...

The charge of this atom type. In the case of oxygen it is 8.

Number of unique atoms
----------------------

::

    DIRAC
    Water
    aug-cc-pVDZ basis (note that a tight p has NOT been added)
    C   2    2  X  Y
            8.    1
    O      .0000000000        0.0000000000        -.2249058930
    ...

Number of symmetry-distinct atoms of this type, which is 1 for oxygen in
water.

Nuclear exponent
----------------

::

    DIRAC
    Water
    aug-cc-pVDZ basis (note that a tight p has NOT been added)
    C   2    2  X  Y
            8.    1                    
    O      .0000000000        0.0000000000        -.2249058930
    ...

Gaussian exponent for the nuclear charge distribution. If blank, then
the default value is used.

Atom name
---------

::

    DIRAC
    Water
    aug-cc-pVDZ basis (note that a tight p has NOT been added)
    C   2    2  X  Y
            8.    1
    O      .0000000000        0.0000000000        -.2249058930
    ...

Is is recommended to use a different name for each atom of the same type
in order to be able to identify them in the output. Usually one labels
symmetry-distinct atoms of the same proton number with indexes. For
example, hydrogen atoms in the water molecule (C1 symmetry):

::

     ...
            1.     2
    H1    1.4523499293        0.0000000000        0.8996235720
    H2   -1.4523499293        0.0000000000        0.8996235720
     ...

Cartesian coordinates
---------------------

::

    DIRAC
    Water
    aug-cc-pVDZ basis (note that a tight p has NOT been added)
    C   2    2  X  Y
            8.    1
    O      .0000000000        0.0000000000        -.2249058930
    ...

The *x*-, *y*-, and *z* coordinates are in Bohr.

::

    DIRAC
    Water
    aug-cc-pVDZ basis (note that a tight p has NOT been added)
    C   2    2  X  Y    
    O      .0000000000        0.0000000000        -.2249058930
    ...

If not blank (usually one places *A* character there) then the molecular
coordinates are read in ångströms instead of Bohr.

The Cartesian coordinates may be given in free format. However, the name
of the atom must still be left four places, and no coordinates must
enter the four first positions.

