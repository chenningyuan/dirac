BSSE calculations with DIRAC
============================

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

Contents
--------

-  `1 Definition of ghost atoms <#Definition_of_ghost_atoms>`_
-  `2 DFT grid <#DFT_grid>`_
-  `3 Symmetry recognition <#Symmetry_recognition>`_

 The basis set superposition error (BSSE) can be estimated by performing
subsystem calculations using ghost atoms which carry basis sets of the
full system.

This tutorial describes how to use ghost atoms in DIRAC and how to avoid
certain pitfalls.

Definition of ghost atoms
-------------------------

In order to replace a real atom by its ghost atom, it is sufficient to
replace the charge by zero charge. It is recommended to also change the
name of the atom. Note that you have to give the ghost atom basis
explicitly. For instance, the input for a ghost atom with a contracted
oxygen aug-cc-pVDZ basis could be given by the following molecule file:

::

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

DFT grid
--------

When you replace an atom by a ghost atom the DFT grid will change. In
order to use the same grid for both the full system and the subsystem
with ghost atoms, you can export the DFT grid of the full system:

::

    *DFT
    .EXPORT
    grid_file

and import this grid back into the calculation with ghost atoms:

::

    *DFT
    .IMPORT
    grid_file

Make sure to copy grid\_file from/to the working directory by using

::

    pam ... -get grid_file ...

and

::

    pam ... -copy grid_file ...

Symmetry recognition
--------------------

If you leave the symmetry text in the molecule specification blank,
DIRAC will try to evaluate the symmetry of the system automatically.
Note that the symmetry recognition uses the atom types and coordinates
but not the respective basis sets. This means that if you use ghost
atoms with different basis sets, you should give the symmetry explicitly
or check the symmetry recognized by DIRAC.

Retrieved from
"`http://wiki.chem.vu.nl/dirac/index.php/BSSE\_calculations\_with\_DIRAC <http://wiki.chem.vu.nl/dirac/index.php/BSSE_calculations_with_DIRAC>`_"

Views
'''''

-  `Article </dirac/index.php/BSSE_calculations_with_DIRAC>`_
-  `Discussion </dirac/index.php?title=Talk:BSSE_calculations_with_DIRAC&action=edit>`_
-  `View
   source </dirac/index.php?title=BSSE_calculations_with_DIRAC&action=edit>`_
-  `History </dirac/index.php?title=BSSE_calculations_with_DIRAC&action=history>`_

Personal tools
''''''''''''''

-  `Log in / create
   account </dirac/index.php?title=Special:Userlogin&returnto=BSSE_calculations_with_DIRAC>`_

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
   here </dirac/index.php/Special:Whatlinkshere/BSSE_calculations_with_DIRAC>`_
-  `Related
   changes </dirac/index.php/Special:Recentchangeslinked/BSSE_calculations_with_DIRAC>`_
-  `Upload file </dirac/index.php/Special:Upload>`_
-  `Special pages </dirac/index.php/Special:Specialpages>`_
-  `Printable
   version </dirac/index.php?title=BSSE_calculations_with_DIRAC&printable=yes>`_
-  `Permanent
   link </dirac/index.php?title=BSSE_calculations_with_DIRAC&oldid=24601>`_

`|Powered by MediaWiki| <http://www.mediawiki.org/>`_

-  This page was last modified 19:59, 25 November 2010.
-  This page has been accessed 1,140 times.
-  `Privacy policy </dirac/index.php/Diracwiki:Privacy_policy>`_
-  `About Diracwiki </dirac/index.php/Diracwiki:About>`_
-  `Disclaimers </dirac/index.php/Diracwiki:General_disclaimer>`_

.. |Powered by
MediaWiki| image:: /dirac/skins/common/images/poweredby_mediawiki_88x31.png
