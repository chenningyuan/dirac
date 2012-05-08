Manual:Basis:General
====================

General information
===================

Choosing input type
-------------------

Presently Dirac can handle two types of files with molecular
coordinates:

-  `.MOL </dirac/index.php/Manual:Basis:MOL>`_
-  `.XYZ </dirac/index.php/Manual:Basis:XYZ>`_

Which file is present, is determined on the basis of whether a dot is
present in the name of the submitted molecule file. If no dot is
present, it is assumed to be a .MOL file, if the name ends with .MOL or
.XYZ dirac will detect this and handle the molecule file accordingly. If
an other extension is detected Dirac will exit (as more filetypes may be
included later). In a .MOL file all the information of the molecule,
like basis set, coordinates, charges etc. is present in the molecule
file. In a .XYZ file only the coordinates of the atoms are given. The
setup of the calculation is given in the .INP file.

Choosing basis sets
-------------------

Apart from the practical aspects involved in setting up a MOLECULE file
for Dirac, one must keep in mind also the more fundamental issue of
choosing a suitable basis set for a given calculation.

The development of basis sets suitable to use in relativistic
calculations reflects the relatively lateness of the field's
development. Because the more consistent efforts in method development
started at about the mid 1980's, it wasn't until well into the late
1990's that the pioneering works of the early and mid 1990's were
substantially complemented and improved upon.

This situations has dramatically improved in recent years, notably with
the work of K. G. Dyall, and Dirac users are strongly advised to use
Dyall's basis sets whenever they are available. These sets follow
roughly the "correlation-consisten" philosophy introduced by Dunning and
coworkers
[1] <#_note-dunninc:cc:1> 
, so they already
contain polarization functions.

At the time of writing, these are avaiable from Ga onwards, and a
detailed description can be found at the dirac website, under the "Basis
Set Repository"
link`[1] <http://dirac.chem.sdu.dk/basisarchives/dyall/index.html>`_.
For convenience, Dirac collects all of Dyall's available basis sets
(those published and areadly in distributable form, as the ones listed
in the repository below, as well as those in preparation) in basis set
files. One must keep in mind that currently these files support only the
use of fully uncontracted basis sets.

The following files are present in Dirac's distribution:

dyall.v2z contains the double-zeta basis dyall.v3z contains the
triple-zeta basis dyall.v4z contains the quadruple-zeta basis dyall.av2z
contains the double-zeta basis (as in dyall.v2z) plus diffuse functions
(p-block only) dyall.av3z contains the triple-zeta basis (as in
dyall.v3z) plus diffuse functions (p-block only) dyall.av4z contains the
quadruple-zeta basis (as in dyall.v4z) plus diffuse functions (p-block
only)

For light elements (say up to Argon) it is advisable to continue using
standard non-relativistic basis set, such as the correlation-consistent
sets of Dunning and coworkers. It is advisable that, in order to have a
balanced description when light and heavy elements are present, that one
uses either contracted or uncontracted sets thoughout.

Apart from Dyall's sets, one can choose several different basis sets
based upon geometric progressions of exponents. One such set is that of
K. Faegri, also available in the basis set library, but with the
drawback that the user needs to extended it by adding polarization
functions.

The procedure used in Dirac is as follows: first, the basis set for the
large component will be read in from the the MOLECULE file; then the
small component basis set is generated automatically via the kinetic
balance prescription, but can (by experienced users) also be specified
explicitly. Try this only if you have studied the theory and know what
the criteria for small component basis sets are, there are many pitfalls
that one needs to be aware of!

<#_ref-dunninc:cc:1_0>`_ T. H. Dunning, *bla*, J. Chem. Phys.
**106**, *0* (1989) `electronic
version <http://dx.doi.org/10.1007/s002140100263>`_.


Manual:Basis:Atomic coordinates
===============================

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


Manual:Basis:Fitting basis set
==============================

**This page describes unreleased functionality. The keywords may not be
available in your version of DIRAC.**

Fitting basis set
=================

In the development version of DIRAC it is now possible to make use of
density fitting. In such calculations, it is mandatory to include a
basis set used to expand the density.

This is done adding a line with the FTSET keyword to the MOLECULE file,
after the coordinates for a given atom type

::

    FTSET [options]

Basically all the options available in the definition of large or small
basis sets are applicable to the fit set (i.e. the use of a set from the
library, explicitly typed basis, even/well-tempered etc). Here are some
examples:

Basis sets from the basis set library
-------------------------------------

::

    DIRAC
    Water
    aug-cc-pVDZ basis (note that a tight p has NOT been added)
    C   2    2  X  Y
            8.    1
    O      .0000000000        0.0000000000        -.2249058930
    LARGE BASIS cc-pVDZ
    FTSET BASIS Ahlrichs-Coulomb-Fit


Manual:Basis:Large component basis set
======================================

Large component basis set
=========================

The specification of the large component basis set follows the atomic
coordinates. There are various possibilities for giving the basis set.

Basis sets from the basis set library
-------------------------------------

::

    DIRAC
    Water
    aug-cc-pVDZ basis (note that a tight p has NOT been added)
    C   2    2  X  Y
            8.    1
    O      .0000000000        0.0000000000        -.2249058930
    LARGE BASIS aug-cc-pVDZ

If you use a aug-cc-p... basis set, you can prefix it with d-, t-, q-,
and it will add the diffuse functions in an even-tempered fashion.

Multiple basis sets from the basis set library
----------------------------------------------

As an alternative to the BASIS option described above, it is possible to
use specify different basis set files via the MULTIBASIS keyword. The
MULTIBASIS option is present to allow easier inclusion of different sets
of diffuse and/or polarization functions to a reference basis set. The
syntax for this keyword is very similar to that of the BASIS:

::

    DIRAC
    Water
    aug-cc-pVDZ basis (note that a tight p has NOT been added)
    C   2    2  X  Y
            8.    1
    O      .0000000000        0.0000000000        -.2249058930
    LARGE MULTIBASIS 2 cc-pVDZ cc-pVDZ-diffuse

After the keyword an integer with the number of files to be read should
be specified. It is followed by the name of the different basis set
files, each separated by a whitespace. The only limitation for the
number of basis set files is that the total lenght of this line should
not exceed Fortran's maximum allowed line size (80 characters).

Explicitly typed basis sets
---------------------------

Explicitly typed basis sets are best described using an explicit
example:

::

    DIRAC
    Water
    aug-cc-pVDZ basis (note that a tight p has NOT been added)
    C   2    2  X  Y
            8.    1
    O      .0000000000        0.0000000000        -.2249058930
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

Highest angular quantum number *l* plus one. In this case it is 3, since
we are using a *spd* basis set.

Number of blocks for each *l*-value. The memory requirements grow
rapidly with the number of basis functions in a block (note for instance
that four *g* functions actually are 60 basis functions, as there are 15
cartesian components of each *g* function). Memory requirements can
therefore be reduced by splitting basis functions of the quantum number
into different blocks. This will, however, decrease the performance of
the integral calculation.

Lines starting with either !, $, or # are comments.

Number of primitive Gaussians in this block.

Number of contracted Gaussians in this block. If a zero is given, an
uncontracted basis set will be assumed, and only orbital exponents need
to be given.

::

    DIRAC
    Water
    aug-cc-pVDZ basis (note that a tight p has NOT been added)
    C   2    2  X  Y
            8.    1
    O      .0000000000        0.0000000000        -.2249058930
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

A single character describing the input format of the basis set in this
block.

The default format is 8F10.4 which will be used if left blank. **Be very
careful when using this default format as it will miss any exponential
parameter standing to the right of the 10 characters**. In this format
the first column is the orbital exponent and the seven last columns are
contraction coefficients. If no numbers are given, a zero is assumed. If
more than 7 contracted functions occur in a given block, the contraction
coefficients may be continued on the next line, but the first column
(where the orbital exponents are given) must then be left blank.

An F or f in the first position (like in the example above) will
indicate that the input is in free format. This will of course require
that all contraction coefficients need to be typed in, as all numbers
need to be present on each line. However, note that this options is
particularly handy together with completely decontracted basis sets, as
described below. Note that the program reads the free format input from
an internal file that is 80 characters long, and no line should
therefore exceed 80 characters.

One may also give the format H or h. This corresponds to high precision
format 4F20.8, where the first column again is reserved for the orbital
exponents, and the three next columns are designated to the contraction
coefficients. If no number is given, a zero is assumed. If there are
more than three contracted orbitals in a given block, the contraction
coefficients may be continued on the next line, though keeping the
column of the orbital exponents blank.

::

    DIRAC
    Water
    aug-cc-pVDZ basis (note that a tight p has NOT been added)
    C   2    2  X  Y
            8.    1
    O      .0000000000        0.0000000000        -.2249058930
    LARGE EXPLICIT  3    1    1    1
    # s functions
    f  10    4    0
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
    f   5    3    0
         17.7000000  0.04301800  0.00000000  0.00000000
          3.8540000  0.22891300  0.00000000  0.00000000
          1.0460000  0.50872800  0.00000000  0.00000000
          0.2753000  0.46053100  1.00000000  0.00000000
          0.0685600  0.00000000  0.00000000  1.00000000
    # d functions
    f   2    2    0
          1.1850000  1.00000000  0.00000000
          0.3320000  0.00000000  1.00000000

Specification of how to generate small component functions by kinetic
balance. If no number or a 0 is given, the small component functions are
generated both upwards and downwards, e.g 
*p* \\rightarrow| *s, d*. 
If
the number is 1, the small component functions are generated upwards,
e.g *p* -->  *d*. 
If the number is 2, the small component functions
are generated downwards, e.g *p* -->  *s*.

For other values no small components functions generated.

The routine that generates the small component basis set by kinetic
balance will delete duplicate functions.

MOLFDIR-type basis sets
-----------------------

::

    DIRAC
    Water
    MOLFDIR basis
    C   2    2  X  Y
            8.    1
    O      .0000000000        0.0000000000        -.2249058930
    LARGE MOLFBAS Oxygen-xyz.bas

The MOLFDIR basis set file(s) (here: Oxygen-xyz.bas) must be copied to
the scratch area, for example using the pam script:

::

    pam ... -copy "Oxygen-xyz.bas ..."

Even-tempered basis sets (geometric progressions)
-------------------------------------------------

The exponents η are generated in an even-tempered series:

η\ :sub:`*N* − *k* + 1`\  = α β\ :sup:`*k* − 1`\ ; *k* = 1, ..., *N*

::

    DIRAC
    Water
    even-tempered basis set
    C   2    2  X  Y
            8.    1
    O      .0000000000        0.0000000000        -.2249058930
    LARGE EVENTEMP 0.05 3.0 11 3 2 1 1
    1..5
    6..11
    7..11
    9..10

::

    DIRAC
    Water
    even-tempered basis set
    C   2    2  X  Y
            8.    1
    O      .0000000000        0.0000000000        -.2249058930
    LARGE EVENTEMP 0.05 3.0 11 3 2 1 1
    1..5
    6..11
    7..11
    9..10

Highest angular quantum number *l* plus one. In this case it is 3, since
we are using a *spd* basis set.

Number of blocks for each *l*-value. The memory requirements grow
rapidly with the number of basis functions in a block (note for instance
that four *g* functions actually are 60 basis functions, as there are 15
cartesian components of each *g* function). Memory requirements can
therefore be reduced by splitting basis functions of the quantum number
into different blocks. This will, however, decrease the performance of
the integral calculation.

::

    DIRAC
    Water
    even-tempered basis set
    C   2    2  X  Y
            8.    1
    O      .0000000000        0.0000000000        -.2249058930
    LARGE EVENTEMP 0.05 3.0 11 3 2 1 1
    1..5
    6..11
    7..11
    9..10

For each block the range of exponent used for that particular block is
given. In this case the first block (*s* functions) will consist of
exponents 1-5, the second block (*s* functions) of exponents 6-11, the
third block (*p* functions) of exponents 7-11, and the last block (*d*
functions) of exponents 9-10.

Well-tempered basis sets
------------------------

The exponents η are generated in an well-tempered series:

η\ :sub:`*N*`\  = α

η\ :sub:`*N* − *k* + 1`\  = η\ :sub:`*N* − *k* + 2`\  
β 
\\bigl[ 1 + | γ
\\bigl( \\frac{k}{N} \\bigr)| \ :sup:`δ`\ 
\\bigr]|; *k* = 1, ..., *N*

::

    DIRAC
    Water
    well-tempered basis set
    C   2    2  X  Y
            8.    1
    O      .0000000000        0.0000000000        -.2249058930
    LARGE WELLTEMP 0.05 2.5 2.0 6.0 11 3 2 1 1
    1..5
    6..11
    7..11
    9..10

For the specification of blocks and exponent ranges see `LARGE
EVENTEMP </dirac/index.php/Manual:Basis:Large_component_basis_set#Even-tempered_basis_sets_.28geometric_progressions.29>`_.

Family basis sets
-----------------

Input for basis sets where the same set of exponents are used for all
functions. This is analogous to the well- and even-tempered basis sets
except that the exponents are not calculated from a formula, but must be
given in the file. These exponents may come from a basis set
optimization with GRASP.

::

    DIRAC
    Water
    family basis set
    C   2    2  X  Y
            8.    1
    O      .0000000000        0.0000000000        -.2249058930
    LARGE FAMILY   10    3    1    1    1
    # exponents
     6665.0000000
     1000.0000000
      228.0000000
       64.7100000
       21.0600000
        7.4950000
        2.7970000
        0.5215000
        0.1596000
        0.0469000
    # ranges
    1..10
    6..10
    8..9

For the specification of blocks and exponent ranges see `LARGE
EVENTEMP </dirac/index.php/Manual:Basis:Large_component_basis_set#Even-tempered_basis_sets_.28geometric_progressions.29>`_.

Dual family basis sets
----------------------

A basis set analogous to the `family basis
set </dirac/index.php/Manual:Basis:Large_component_basis_set#Family_basis_sets>`_,
except one set of exponents are used for *s*, *d*, *g*, ... functions
and another set is used for *p*, *f*, *h*, ... functions.

::

    DIRAC
    Water
    dual family basis set
    C   2    2  X  Y
            8.    1
    O      .0000000000        0.0000000000        -.2249058930
    LARGE DUALFAMILY   10    5    3    1    1    1
    # s, d, g, ...
     6665.0000000
     1000.0000000
      228.0000000
       64.7100000
       21.0600000
        7.4950000
        2.7970000
        0.5215000
        0.1596000
        0.0469000
    # p, f, h, ...
        9.4390000
        2.0020000
        0.5456000
        0.1517000
        0.0404100
    # ranges
    1..10
    1..5
    8..9

For the specification of blocks and exponent ranges see 
`LARGE
EVENTEMP </dirac/index.php/Manual:Basis:Large_component_basis_set#Even-tempered_basis_sets_.28geometric_progressions.29>`_.


Manual:Basis:MOL
================

DIRAC Molecule File Structure
-----------------------------

In this section we provide more information about the MOLECULE format as
is used in Dirac. This file is divided in sections, with the following
being mandatory:

-  `Title and
   symmetry </dirac/index.php/Manual:Basis:Title_and_symmetry>`_
-  `Atomic
   coordinates </dirac/index.php/Manual:Basis:Atomic_coordinates>`_
-  `Large component basis
   set </dirac/index.php/Manual:Basis:Large_component_basis_set>`_

It is possible to define explicitly the small component basis set to be
use with an additional section:

-  `Small component basis
   set </dirac/index.php/Manual:Basis:Small_component_basis_set>`_

but this is however not generally advised, given that the small
component basis set is more conveniently generated via the kinetic
balance procedure.

In the development version(s) of DIRAC it is also possible to define the
so-called fit sets via an additional section:

-  `Fitting basis
   set </dirac/index.php/Manual:Basis:Fitting_basis_set>`_



Small component basis set
=========================

After the large component basis set the small component basis set *can*
be specified. If nothing is specified it is equivalent to specifying
SMALL KINBAL (see below). There are three possibilities for giving the
basis set.

By kinetic balance
------------------

::

    DIRAC
    Water
    aug-cc-pVDZ basis (note that a tight p has NOT been added)
    C   2    2  X  Y
            8.    1
    O      .0000000000        0.0000000000        -.2249058930
    LARGE BASIS aug-cc-pVDZ
    SMALL KINBAL

The small component basis is generated using kinetic balance. The small
component basis set will *always* be uncontracted, independent of
whether the large component basis set is contracted or not.

Explicitly typed basis sets
---------------------------

::

    DIRAC
    Water
    aug-cc-pVDZ basis (note that a tight p has NOT been added)
    C   2    2  X  Y
            8.    1
    O      .0000000000        0.0000000000        -.2249058930
    LARGE BASIS aug-cc-pVDZ
    SMALL EXPLICIT  3    1    1    1
    ...

Syntax is analogous to `LARGE
EXPLICIT </dirac/index.php/Manual:Basis:Large_component_basis_set#Explicitly_typed_basis_sets>`_.

MOLFDIR-type basis sets
-----------------------

::

    DIRAC
    Water
    MOLFDIR basis
    C   2    2  X  Y
            8.    1
    O      .0000000000        0.0000000000        -.2249058930
    LARGE MOLFBAS Oxygen-xyz.bas
    SMALL MOLFBAS Oxygen-xyz.bas

The MOLFDIR basis set file(s) (here: Oxygen-xyz.bas) must be copied to
the scratch area, for example using the pam script:

::

    pam ... --copy "Oxygen-xyz.bas ..."


Title and symmetry
==================

The detailed input is best described using an example:

::

    DIRAC
    Water
    aug-cc-pVDZ basis (note that a tight p has NOT been added)
    C   2    2  X  Y
            8.    1
    O      .0000000000        0.0000000000        -.2249058930
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

First line
----------

We start to focus at the first line in the input:

::

    DIRAC 
    Water
    aug-cc-pVDZ basis (note that a tight p has NOT been added)
    C   2    2  X  Y

A keyword that is not used in DIRAC, it can be any sequence of
characters.

Next two lines
--------------

Next two lines are arbitrary.

::

    DIRAC
    Water                                                                  
    aug-cc-pVDZ basis (note that a tight p has NOT been added)             
    C   2    2  X  Y

We recommend to place there short descriptions of basis sets, system, or
any other usefull information.

Cartesian Gaussians
-------------------

::

    DIRAC
    Water
    aug-cc-pVDZ basis (note that a tight p has NOT been added)
    C   2    2  X  Y

Must be set to C to indicate that Cartesian Gaussians are to be used. If
not, DIRAC will issue a warning and continue assuming Cartesian
Gaussians. Note that by default transformation to spherical harmonic
basis (modified for the small components) is done in the transformation
to the orthogonal basis. Spherical Gaussians (S) are allowed in
2-component calculations.

Number of atomic types
----------------------

Specify the number of atomic types, which is 2 in this example (oxygen
and hydrogen).

::

    DIRAC
    Water
    aug-cc-pVDZ basis (note that a tight p has NOT been added)
    C   2    2  X  Y

Charge of the system
--------------------

::

    DIRAC
    Water
    aug-cc-pVDZ basis (note that a tight p has NOT been added)
    C   2    2  X  Y

Specify the charge of the system. This will be used to determine the
initial Hartree-Fock/Kohn-Sham occupation, if it is not explicitly
specified. A zero (no charge) is assumed if no number is given here.

The charge functionality is supressed by explicitly setting the
`occupation </dirac/index.php/Manual:WAVE_FUNCTION:DHFCAL#Occupation>`_.

Number of symmetry generators
-----------------------------

::

    DIRAC
    Water
    aug-cc-pVDZ basis (note that a tight p has NOT been added)
    C   2    2  X  Y

Specify the number of generators of the symmetry group. For
*C*\ :sub:`2*v*`\  there are two generators, in this example
σ\ :sub:`*x*`\  and σ\ :sub:`*y*`\ . You must specify explicitly a zero
0 if you want to run the calculation without symmetry. Be ware that in
case of symmetry you should only specify a *symmetry-unique* list of
atoms, as the program will use the generators to find the full list of
atoms. In case of linear molecules one can specify Dinfh or Cinfv so
that the full symmetry group D\_{\\infty h} or C\_{\\infty v} is
employed in some of the modules. This symmetry is not implemented in all
modules, so that you may also see output referring to the subgroups
*D*\ :sub:`2*h*`\  or *C*\ :sub:`2*v*`\ , respectively. If the keyword
is left blank and all atoms are specified explicitly, then DIRAC will
analyze the geometry and find the appropriate symmetry group by itself.
In doing so it may rotate the molecule.

Symmetry group string
---------------------

::

    DIRAC
    Water
    aug-cc-pVDZ basis (note that a tight p has NOT been added)
    C   2    2  X  Y   

In general, the symmetry handling will be restricted to the binary
groups, that is *D*\ :sub:`2*h*`\  and subgroups, which means that a
symmetry operation acting on the main axes (x,y,z) will at most reverse
their direction. A group generator is therefore identified by a
3-character string that specifies the axes reversed under its operation.
In this example the X and Y are given to indicate the generators
σ\ :sub:`*x*`\  and σ\ :sub:`*y*`\ .

Used DIRAC binary groups
------------------------

Examples for the eight binary groups:

**Group**

**Operations**

*D*\ :sub:`2*h*`\ 

::

    C   ?    3  Z  Y  X

σ(*x**y*),σ(*x**z*),σ(*y**z*)

*D*\ :sub:`2`\ 

::

    C   ?    2XY  YZ   

*C*\ :sub:`2`\ (*z*),*C*\ :sub:`2`\ (*x*)

*C*\ :sub:`2*v*`\ 

::

    C   ?    2 Y X     

σ(*x**z*),σ(*y**z*)

*C*\ :sub:`2*h*`\ 

::

    C   ?    2  ZXYZ   

σ(*x**y*),*i*

*C*\ :sub:`2`\ 

::

    C   ?    1XY       

*C*\ :sub:`2`\ (*z*)

*C*\ :sub:`*s*`\ 

::

    C   ?    1  Z      

σ(*x**y*)

*C*\ :sub:`*i*`\ 

::

    C   ?    1XYZ      

*i*

*C*\ :sub:`1`\ 

::

    C   ?    0         

Units of coordinates
--------------------

::

    DIRAC
    Water
    aug-cc-pVDZ basis (note that a tight p has NOT been added)
    C   2    2  X  Y    

If not blank then the molecular coordinates are read in ångström instead
of bohr.

Neglect integrals threshold
---------------------------

::

    DIRAC
    Water
    aug-cc-pVDZ basis (note that a tight p has NOT been added)
    C   2    2  X  Y              

Threshold for neglect of final integrals. Default is 1.0D-15, which will
be used if no input is given here. A threshold of 1.0D-15 will give
integrals correct to approximately 1.0D-13. Separate thresholds for SL
and SS integrals may be specified in the input file (see
`.THRFAC </dirac/index.php/Manual:INTEGRALS:TWOINT#.THRFAC>`_).


