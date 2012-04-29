Manual:DIRRCI:CIROOT
====================

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

Contents
--------

-  `1 &CIROOT -- Select the state to converge
   on <#.26CIROOT_--_Select_the_state_to_converge_on>`_

   -  `1.1 IREPNA <#IREPNA>`_
   -  `1.2 NROOTS <#NROOTS>`_
   -  `1.3 ISTART <#ISTART>`_
   -  `1.4 NSEL(NROOTS) <#NSEL.28NROOTS.29>`_
   -  `1.5 SELECT <#SELECT>`_

&CIROOT -- Select the state to converge on
==========================================

IREPNA
------

Abelian symmetry group of the desired state(s).

*Default: First abelian symmetry in the list*

NROOTS
------

Number of states to optimize on.

*Default:*

::

    NROOTS=1

`|image0| </dirac/index.php/Image:Dirac.png>`_

**Advanced options**

ISTART
------

Start vector method.

COSCI start vectors:

::

    ISTART=1

Determinant with lowest eigenvalue:

::

    ISTART=2

First (reference) determinant:

::

    ISTART=3

If the COSCI start vectors are not available, then the default will be
2.

NSEL(NROOTS)
------------

Rank number of states to be optimized.

*Default:*

::

    NSEL=1,2,3,...

SELECT
------

Select wave functions on basis of largest overlap with the start wave
function.

*Default:*

::

    SELECT=F

Retrieved from
"`http://wiki.chem.vu.nl/dirac/index.php/Manual:DIRRCI:CIROOT <http://wiki.chem.vu.nl/dirac/index.php/Manual:DIRRCI:CIROOT>`_"

Views
'''''

-  `Article </dirac/index.php/Manual:DIRRCI:CIROOT>`_
-  `Discussion </dirac/index.php?title=Talk:Manual:DIRRCI:CIROOT&action=edit>`_
-  `View
   source </dirac/index.php?title=Manual:DIRRCI:CIROOT&action=edit>`_
-  `History </dirac/index.php?title=Manual:DIRRCI:CIROOT&action=history>`_

Personal tools
''''''''''''''

-  `Log in / create
   account </dirac/index.php?title=Special:Userlogin&returnto=Manual:DIRRCI:CIROOT>`_

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
   here </dirac/index.php/Special:Whatlinkshere/Manual:DIRRCI:CIROOT>`_
-  `Related
   changes </dirac/index.php/Special:Recentchangeslinked/Manual:DIRRCI:CIROOT>`_
-  `Upload file </dirac/index.php/Special:Upload>`_
-  `Special pages </dirac/index.php/Special:Specialpages>`_
-  `Printable
   version </dirac/index.php?title=Manual:DIRRCI:CIROOT&printable=yes>`_
-  `Permanent
   link </dirac/index.php?title=Manual:DIRRCI:CIROOT&oldid=2401>`_

`|Powered by MediaWiki| <http://www.mediawiki.org/>`_

-  This page was last modified 15:15, 5 March 2008.
-  This page has been accessed 909 times.
-  `Privacy policy </dirac/index.php/Diracwiki:Privacy_policy>`_
-  `About Diracwiki </dirac/index.php/Diracwiki:About>`_
-  `Disclaimers </dirac/index.php/Diracwiki:General_disclaimer>`_

.. |image0| image:: /dirac/images/thumb/9/95/Dirac.png/41px-Dirac.png
.. |Powered by
MediaWiki| image:: /dirac/skins/common/images/poweredby_mediawiki_88x31.png
