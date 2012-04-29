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

Retrieved from
"`http://wiki.chem.vu.nl/dirac/index.php/Manual:Basis:Atomic\_coordinates <http://wiki.chem.vu.nl/dirac/index.php/Manual:Basis:Atomic_coordinates>`_"

Views
'''''

-  `Article </dirac/index.php/Manual:Basis:Atomic_coordinates>`_
-  `Discussion </dirac/index.php?title=Talk:Manual:Basis:Atomic_coordinates&action=edit>`_
-  `View
   source </dirac/index.php?title=Manual:Basis:Atomic_coordinates&action=edit>`_
-  `History </dirac/index.php?title=Manual:Basis:Atomic_coordinates&action=history>`_

Personal tools
''''''''''''''

-  `Log in / create
   account </dirac/index.php?title=Special:Userlogin&returnto=Manual:Basis:Atomic_coordinates>`_

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
   here </dirac/index.php/Special:Whatlinkshere/Manual:Basis:Atomic_coordinates>`_
-  `Related
   changes </dirac/index.php/Special:Recentchangeslinked/Manual:Basis:Atomic_coordinates>`_
-  `Upload file </dirac/index.php/Special:Upload>`_
-  `Special pages </dirac/index.php/Special:Specialpages>`_
-  `Printable
   version </dirac/index.php?title=Manual:Basis:Atomic_coordinates&printable=yes>`_
-  `Permanent
   link </dirac/index.php?title=Manual:Basis:Atomic_coordinates&oldid=14311>`_

`|Powered by MediaWiki| <http://www.mediawiki.org/>`_

-  This page was last modified 19:12, 26 October 2008.
-  This page has been accessed 1,951 times.
-  `Privacy policy </dirac/index.php/Diracwiki:Privacy_policy>`_
-  `About Diracwiki </dirac/index.php/Diracwiki:About>`_
-  `Disclaimers </dirac/index.php/Diracwiki:General_disclaimer>`_

.. |Powered by
MediaWiki| image:: /dirac/skins/common/images/poweredby_mediawiki_88x31.png
