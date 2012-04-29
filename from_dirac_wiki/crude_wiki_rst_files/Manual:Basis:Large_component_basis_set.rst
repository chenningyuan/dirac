Manual:Basis:Large component basis set
======================================

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

Contents
--------

-  `1 Large component basis set <#Large_component_basis_set>`_

   -  `1.1 Basis sets from the basis set
      library <#Basis_sets_from_the_basis_set_library>`_
   -  `1.2 Multiple basis sets from the basis set
      library <#Multiple_basis_sets_from_the_basis_set_library>`_
   -  `1.3 Explicitly typed basis sets <#Explicitly_typed_basis_sets>`_
   -  `1.4 MOLFDIR-type basis sets <#MOLFDIR-type_basis_sets>`_
   -  `1.5 Even-tempered basis sets (geometric
      progressions) <#Even-tempered_basis_sets_.28geometric_progressions.29>`_
   -  `1.6 Well-tempered basis sets <#Well-tempered_basis_sets>`_
   -  `1.7 Family basis sets <#Family_basis_sets>`_
   -  `1.8 Dual family basis sets <#Dual_family_basis_sets>`_

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
generated both upwards and downwards, e.g *p* |\\rightarrow| *s, d*. If
the number is 1, the small component functions are generated upwards,
e.g *p* |image1| *d*. If the number is 2, the small component functions
are generated downwards, e.g *p* |image2| *s*.

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

η\ :sub:`*N* − *k* + 1`\  = η\ :sub:`*N* − *k* + 2`\  β |\\bigl[ 1 + | γ
|\\bigl( \\frac{k}{N} \\bigr)| \ :sup:`δ`\  |\\bigr]|; *k* = 1, ..., *N*

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

For the specification of blocks and exponent ranges see `LARGE
EVENTEMP </dirac/index.php/Manual:Basis:Large_component_basis_set#Even-tempered_basis_sets_.28geometric_progressions.29>`_.

