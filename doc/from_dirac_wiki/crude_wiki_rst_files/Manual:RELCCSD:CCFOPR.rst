Manual:RELCCSD:CCFOPR
=====================

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

Contents
--------

-  `1 &CCFOPR -- First-order
   properties <#.26CCFOPR_--_First-order_properties>`_

   -  `1.1 DOMP2G <#DOMP2G>`_
   -  `1.2 DOMP2GOLD <#DOMP2GOLD>`_
   -  `1.3 NEOPER <#NEOPER>`_
   -  `1.4 NAMEE <#NAMEE>`_
   -  `1.5 NTOL <#NTOL>`_
   -  `1.6 MAXIT <#MAXIT>`_
   -  `1.7 MAXDIM <#MAXDIM>`_

&CCFOPR -- First-order properties
=================================

Only MP2 first order properties are currently available. They are
computed by forming a relaxed density matrix, which is done (partly) in
AO basis to obtain all relaxation contributions. The properties that are
evaluated are those specified in the expectation value part of the
input.

DOMP2G
------

Calculate the MP2 gradient, that is, expectation values using the
relaxed one-electron density matrix. Note that this is not the molecular
gradient, but other first order properties as specified in the
expectation value part of the input.

*Default:*

::

    DOMP2G=F

`|image0| </dirac/index.php/Image:Dirac.png>`_

**Advanced options for the old formalism**

DOMP2GOLD
---------

Calculate the MP2 gradient via the old (approximate) MO formalism.

*Default:*

::

    DOMP2GOLD=F

NEOPER
------

Number of first order property operators |\\langle E \\rangle| (only
relevant for the old formalism).

*Default:*

::

    NEOPER=1

NAMEE
-----

Names of *E* operators in |image2|.

*Arguments:* Character (A8) (NAMEE(I),I = 1,NEOPER).

*Default:*

::

    NAMEE='ZDIPLEN'

NTOL
----

Convergence tolerance of the DIIS error vector (10\ :sup:`--NTOL`\ ).

*Default:*

::

    NTOL=12

MAXIT
-----

Maximum number of iterations to solve the Z-vector equations.

*Default:*

::

    MAXIT=30

MAXDIM
------

Maximum dimension of the reduced space.

*Default:*

::

    MAXDIM=8

Retrieved from
"`http://wiki.chem.vu.nl/dirac/index.php/Manual:RELCCSD:CCFOPR <http://wiki.chem.vu.nl/dirac/index.php/Manual:RELCCSD:CCFOPR>`_"

Views
'''''

-  `Article </dirac/index.php/Manual:RELCCSD:CCFOPR>`_
-  `Discussion </dirac/index.php?title=Talk:Manual:RELCCSD:CCFOPR&action=edit>`_
-  `View
   source </dirac/index.php?title=Manual:RELCCSD:CCFOPR&action=edit>`_
-  `History </dirac/index.php?title=Manual:RELCCSD:CCFOPR&action=history>`_

Personal tools
''''''''''''''

-  `Log in / create
   account </dirac/index.php?title=Special:Userlogin&returnto=Manual:RELCCSD:CCFOPR>`_

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
   here </dirac/index.php/Special:Whatlinkshere/Manual:RELCCSD:CCFOPR>`_
-  `Related
   changes </dirac/index.php/Special:Recentchangeslinked/Manual:RELCCSD:CCFOPR>`_
-  `Upload file </dirac/index.php/Special:Upload>`_
-  `Special pages </dirac/index.php/Special:Specialpages>`_
-  `Printable
   version </dirac/index.php?title=Manual:RELCCSD:CCFOPR&printable=yes>`_
-  `Permanent
   link </dirac/index.php?title=Manual:RELCCSD:CCFOPR&oldid=2621>`_

`|Powered by MediaWiki| <http://www.mediawiki.org/>`_

-  This page was last modified 09:58, 6 March 2008.
-  This page has been accessed 756 times.
-  `Privacy policy </dirac/index.php/Diracwiki:Privacy_policy>`_
-  `About Diracwiki </dirac/index.php/Diracwiki:About>`_
-  `Disclaimers </dirac/index.php/Diracwiki:General_disclaimer>`_

.. |image0| image:: /dirac/images/thumb/9/95/Dirac.png/41px-Dirac.png
.. |\\langle E
\\rangle| image:: /dirac/images/math/c/3/0/c30634d1784b865a08c9f74f65bfb369.png
.. |image2| image:: /dirac/images/math/c/3/0/c30634d1784b865a08c9f74f65bfb369.png
.. |Powered by
MediaWiki| image:: /dirac/skins/common/images/poweredby_mediawiki_88x31.png
