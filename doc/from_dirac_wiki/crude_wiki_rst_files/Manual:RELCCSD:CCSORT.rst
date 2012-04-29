Manual:RELCCSD:CCSORT
=====================

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

Contents
--------

-  `1 &CCSORT -- Sorting module <#.26CCSORT_--_Sorting_module>`_

   -  `1.1 USEOE <#USEOE>`_
   -  `1.2 NORECMP <#NORECMP>`_

&CCSORT -- Sorting module
=========================

USEOE
-----

Use orbital energies from the HF calculation instead of the diagonal
elements of the reconstructed Fock matrix. Putting USEOE=F may solve
problems with low-symmetry open shell systems, in which a Kamers pair is
split over the occupied and virtual space, giving an orbital energy
difference of zero in the default definition.

*Default:*

::

    USEOE=T

NORECMP
-------

Do not recompute the Fock matrix for the determinant selected as
reference for the coupled cluster expansion. This is usually done as a
check on the numerical accuracy of the MO transformation. If you specify
both USEOE and NORECMP the reference Hamiltonian is given by the
original orbital energies (taken from the SCF program) and the
transformed 2-electron integrals. This makes it possible to apply the
molecular mean field approximation.

*Default:*

::

    NORECMP=F

Retrieved from
"`http://wiki.chem.vu.nl/dirac/index.php/Manual:RELCCSD:CCSORT <http://wiki.chem.vu.nl/dirac/index.php/Manual:RELCCSD:CCSORT>`_"

Views
'''''

-  `Article </dirac/index.php/Manual:RELCCSD:CCSORT>`_
-  `Discussion </dirac/index.php?title=Talk:Manual:RELCCSD:CCSORT&action=edit>`_
-  `View
   source </dirac/index.php?title=Manual:RELCCSD:CCSORT&action=edit>`_
-  `History </dirac/index.php?title=Manual:RELCCSD:CCSORT&action=history>`_

Personal tools
''''''''''''''

-  `Log in / create
   account </dirac/index.php?title=Special:Userlogin&returnto=Manual:RELCCSD:CCSORT>`_

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
   here </dirac/index.php/Special:Whatlinkshere/Manual:RELCCSD:CCSORT>`_
-  `Related
   changes </dirac/index.php/Special:Recentchangeslinked/Manual:RELCCSD:CCSORT>`_
-  `Upload file </dirac/index.php/Special:Upload>`_
-  `Special pages </dirac/index.php/Special:Specialpages>`_
-  `Printable
   version </dirac/index.php?title=Manual:RELCCSD:CCSORT&printable=yes>`_
-  `Permanent
   link </dirac/index.php?title=Manual:RELCCSD:CCSORT&oldid=29011>`_

`|Powered by MediaWiki| <http://www.mediawiki.org/>`_

-  This page was last modified 15:26, 18 May 2011.
-  This page has been accessed 918 times.
-  `Privacy policy </dirac/index.php/Diracwiki:Privacy_policy>`_
-  `About Diracwiki </dirac/index.php/Diracwiki:About>`_
-  `Disclaimers </dirac/index.php/Diracwiki:General_disclaimer>`_

.. |Powered by
MediaWiki| image:: /dirac/skins/common/images/poweredby_mediawiki_88x31.png
