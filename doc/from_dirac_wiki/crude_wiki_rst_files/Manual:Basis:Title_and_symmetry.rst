Manual:Basis:Title and symmetry
===============================

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

Contents
--------

-  `1 Title and symmetry <#Title_and_symmetry>`_

   -  `1.1 First line <#First_line>`_
   -  `1.2 Next two lines <#Next_two_lines>`_
   -  `1.3 Cartesian Gaussians <#Cartesian_Gaussians>`_
   -  `1.4 Number of atomic types <#Number_of_atomic_types>`_
   -  `1.5 Charge of the system <#Charge_of_the_system>`_
   -  `1.6 Number of symmetry
      generators <#Number_of_symmetry_generators>`_
   -  `1.7 Symmetry group string <#Symmetry_group_string>`_
   -  `1.8 Used DIRAC binary groups <#Used_DIRAC_binary_groups>`_
   -  `1.9 Units of coordinates <#Units_of_coordinates>`_
   -  `1.10 Neglect integrals threshold <#Neglect_integrals_threshold>`_

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
that the full symmetry group |D\_{\\infty h}| or |C\_{\\infty v}| is
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

Retrieved from
"`http://wiki.chem.vu.nl/dirac/index.php/Manual:Basis:Title\_and\_symmetry <http://wiki.chem.vu.nl/dirac/index.php/Manual:Basis:Title_and_symmetry>`_"

Views
'''''

-  `Article </dirac/index.php/Manual:Basis:Title_and_symmetry>`_
-  `Discussion </dirac/index.php?title=Talk:Manual:Basis:Title_and_symmetry&action=edit>`_
-  `View
   source </dirac/index.php?title=Manual:Basis:Title_and_symmetry&action=edit>`_
-  `History </dirac/index.php?title=Manual:Basis:Title_and_symmetry&action=history>`_

Personal tools
''''''''''''''

-  `Log in / create
   account </dirac/index.php?title=Special:Userlogin&returnto=Manual:Basis:Title_and_symmetry>`_

` </dirac/index.php/Dirac_Program>`_

Navigation
''''''''''

-  `Home </dirac/index.php/Dirac_Program>`_
-  `Documentation </dirac/index.php/Documentation>`_
-  `Reference manual </dirac/index.php/Manual>`_
-  `Discussion forum </dirac/index.php/DiscussionForum>`_
-  `News Archive </dirac/index.php/News_archive>`_
-  `Basis sets </dirac/index.php/Basis_Sets>`_

Search
''''''

Toolbox
'''''''

-  `What links
   here </dirac/index.php/Special:Whatlinkshere/Manual:Basis:Title_and_symmetry>`_
-  `Related
   changes </dirac/index.php/Special:Recentchangeslinked/Manual:Basis:Title_and_symmetry>`_
-  `Upload file </dirac/index.php/Special:Upload>`_
-  `Special pages </dirac/index.php/Special:Specialpages>`_
-  `Printable
   version </dirac/index.php?title=Manual:Basis:Title_and_symmetry&printable=yes>`_
-  `Permanent
   link </dirac/index.php?title=Manual:Basis:Title_and_symmetry&oldid=23451>`_

`|Powered by MediaWiki| <http://www.mediawiki.org/>`_

-  This page was last modified 14:14, 7 October 2010.
-  This page has been accessed 4,011 times.
-  `Privacy policy </dirac/index.php/Diracwiki:Privacy_policy>`_
-  `About Diracwiki </dirac/index.php/Diracwiki:About>`_
-  `Disclaimers </dirac/index.php/Diracwiki:General_disclaimer>`_

.. |D\_{\\infty
h}| image:: /dirac/images/math/1/2/f/12fb65a005c771345f0bb591f3e8e9bd.png
.. |C\_{\\infty
v}| image:: /dirac/images/math/c/d/4/cd4658c2ba79b9b150e5990a51b5cea1.png
.. |Powered by
MediaWiki| image:: /dirac/skins/common/images/poweredby_mediawiki_88x31.png
