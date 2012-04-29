Manual:RELCCSD:CCENER
=====================

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

Contents
--------

-  `1 &CCENER -- Energy module <#.26CCENER_--_Energy_module>`_

   -  `1.1 DOMP2 <#DOMP2>`_
   -  `1.2 DOCCSD <#DOCCSD>`_
   -  `1.3 DOCCSDT <#DOCCSDT>`_
   -  `1.4 NTOL <#NTOL>`_
   -  `1.5 MAXIT <#MAXIT>`_
   -  `1.6 MAXDIM <#MAXDIM>`_
   -  `1.7 NOCCS <#NOCCS>`_
   -  `1.8 NOCCD <#NOCCD>`_

&CCENER -- Energy module
========================

DOMP2
-----

Calculate the MP2 energy.

*Default:*

::

    DOMP2=T

DOCCSD
------

Calculate the CCSD energy.

*Default:*

::

    DOCCSD=T

DOCCSDT
-------

Calculate the CCSD(T) energy. The program calculates three different
non-iterative triple corrections
\ :sup:``[1] <#_note-visscher:relccsd1>`_`\ .

*Default:*

::

    DOCCSDT=T

`|image0| </dirac/index.php/Image:Dirac.png>`_

**Advanced options**

NTOL
----

Convergence tolerance of the DIIS error vector (10\ :sup:`--NTOL`\ ).

*Default:*

::

    NTOL=12

MAXIT
-----

Maximum number of DIIS iterations to solve the CCSD equations.

*Default:*

::

    MAXIT=30

MAXDIM
------

Maximum dimension of the reduced space used in the DIIS extrapolation.

*Default:*

::

    MAXDIM=8

`|image1| </dirac/index.php/Image:HJAaJ.png>`_

**Programmers options**

NOCCS
-----

Ignore single excitations so that the CCD formalism is used. Only for
test purposes.

*Default:*

::

    NOCCS=F

NOCCD
-----

Ignore double excitations so that the CCS formalism is used. Only for
test purposes.

*Default:*

::

    NOCCD=F

Retrieved from
"`http://wiki.chem.vu.nl/dirac/index.php/Manual:RELCCSD:CCENER <http://wiki.chem.vu.nl/dirac/index.php/Manual:RELCCSD:CCENER>`_"

Views
'''''

-  `Article </dirac/index.php/Manual:RELCCSD:CCENER>`_
-  `Discussion </dirac/index.php?title=Talk:Manual:RELCCSD:CCENER&action=edit>`_
-  `View
   source </dirac/index.php?title=Manual:RELCCSD:CCENER&action=edit>`_
-  `History </dirac/index.php?title=Manual:RELCCSD:CCENER&action=history>`_

Personal tools
''''''''''''''

-  `Log in / create
   account </dirac/index.php?title=Special:Userlogin&returnto=Manual:RELCCSD:CCENER>`_

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
   here </dirac/index.php/Special:Whatlinkshere/Manual:RELCCSD:CCENER>`_
-  `Related
   changes </dirac/index.php/Special:Recentchangeslinked/Manual:RELCCSD:CCENER>`_
-  `Upload file </dirac/index.php/Special:Upload>`_
-  `Special pages </dirac/index.php/Special:Specialpages>`_
-  `Printable
   version </dirac/index.php?title=Manual:RELCCSD:CCENER&printable=yes>`_
-  `Permanent
   link </dirac/index.php?title=Manual:RELCCSD:CCENER&oldid=4601>`_

`|Powered by MediaWiki| <http://www.mediawiki.org/>`_

-  This page was last modified 11:36, 12 March 2008.
-  This page has been accessed 1,138 times.
-  `Privacy policy </dirac/index.php/Diracwiki:Privacy_policy>`_
-  `About Diracwiki </dirac/index.php/Diracwiki:About>`_
-  `Disclaimers </dirac/index.php/Diracwiki:General_disclaimer>`_

.. |image0| image:: /dirac/images/thumb/9/95/Dirac.png/41px-Dirac.png
.. |image1| image:: /dirac/images/thumb/6/65/HJAaJ.png/41px-HJAaJ.png
.. |Powered by
MediaWiki| image:: /dirac/skins/common/images/poweredby_mediawiki_88x31.png
