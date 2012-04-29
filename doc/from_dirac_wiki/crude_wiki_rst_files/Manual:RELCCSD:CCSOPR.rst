Manual:RELCCSD:CCSOPR
=====================

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

Contents
--------

-  `1 &CCSOPR -- Second-order
   properties <#.26CCSOPR_--_Second-order_properties>`_

   -  `1.1 NAOPER <#NAOPER>`_
   -  `1.2 NBOPER <#NBOPER>`_
   -  `1.3 NFREQ <#NFREQ>`_
   -  `1.4 NAMEA <#NAMEA>`_
   -  `1.5 NAMEB <#NAMEB>`_
   -  `1.6 EFREQ <#EFREQ>`_
   -  `1.7 NTOL <#NTOL>`_
   -  `1.8 MAXIT <#MAXIT>`_
   -  `1.9 MAXDIM <#MAXDIM>`_

&CCSOPR -- Second-order properties
==================================

Only the SCF wave function can be used at the moment (RPA). The names of
the property integrals should first be defined in the
`.PRPTRA </dirac/index.php/Manual:MOLTRA#.PRPTRA>`_ section of
`\*\*MOLTRA </dirac/index.php/Manual:MOLTRA>`_.

NAOPER
------

Number of A operators in |\\langle\\langle A ; B
\\rangle\\rangle\_\\omega|.

*Default:*

::

    NAOPER=1

NBOPER
------

Number of B operators in |image1|.

*Default:*

::

    NBOPER=1

NFREQ
-----

Number of frequencies in |image2|.

*Default:*

::

    NFREQ=1

NAMEA
-----

Names of A operators in |image3|.

*Arguments:* Character (A8) (NAMEA(I),I = 1,NAOPER).

*Default:*

::

    NAMEA='ZDIPLEN'

NAMEB
-----

Names of B operators in |image4|.

*Arguments:* Character (A8) (NAMEB(I),I = 1,NBOPER).

*Default:*

::

    NAMEB='ZDIPLEN'

EFREQ
-----

Values of frequencies in |image5|.

*Arguments:* Real (EFREQ(I),I = 1,NFREQ).

*Default:*

::

    EFREQ=0.0

`|image6| </dirac/index.php/Image:Dirac.png>`_

**Advanced options**

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
"`http://wiki.chem.vu.nl/dirac/index.php/Manual:RELCCSD:CCSOPR <http://wiki.chem.vu.nl/dirac/index.php/Manual:RELCCSD:CCSOPR>`_"

Views
'''''

-  `Article </dirac/index.php/Manual:RELCCSD:CCSOPR>`_
-  `Discussion </dirac/index.php?title=Talk:Manual:RELCCSD:CCSOPR&action=edit>`_
-  `View
   source </dirac/index.php?title=Manual:RELCCSD:CCSOPR&action=edit>`_
-  `History </dirac/index.php?title=Manual:RELCCSD:CCSOPR&action=history>`_

Personal tools
''''''''''''''

-  `Log in / create
   account </dirac/index.php?title=Special:Userlogin&returnto=Manual:RELCCSD:CCSOPR>`_

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
   here </dirac/index.php/Special:Whatlinkshere/Manual:RELCCSD:CCSOPR>`_
-  `Related
   changes </dirac/index.php/Special:Recentchangeslinked/Manual:RELCCSD:CCSOPR>`_
-  `Upload file </dirac/index.php/Special:Upload>`_
-  `Special pages </dirac/index.php/Special:Specialpages>`_
-  `Printable
   version </dirac/index.php?title=Manual:RELCCSD:CCSOPR&printable=yes>`_
-  `Permanent
   link </dirac/index.php?title=Manual:RELCCSD:CCSOPR&oldid=4501>`_

`|Powered by MediaWiki| <http://www.mediawiki.org/>`_

-  This page was last modified 07:36, 12 March 2008.
-  This page has been accessed 808 times.
-  `Privacy policy </dirac/index.php/Diracwiki:Privacy_policy>`_
-  `About Diracwiki </dirac/index.php/Diracwiki:About>`_
-  `Disclaimers </dirac/index.php/Diracwiki:General_disclaimer>`_

.. |\\langle\\langle A ; B
\\rangle\\rangle\_\\omega| image:: /dirac/images/math/5/d/1/5d1ccb648774ad117a597fb5cb7ad174.png
.. |image1| image:: /dirac/images/math/5/d/1/5d1ccb648774ad117a597fb5cb7ad174.png
.. |image2| image:: /dirac/images/math/5/d/1/5d1ccb648774ad117a597fb5cb7ad174.png
.. |image3| image:: /dirac/images/math/5/d/1/5d1ccb648774ad117a597fb5cb7ad174.png
.. |image4| image:: /dirac/images/math/5/d/1/5d1ccb648774ad117a597fb5cb7ad174.png
.. |image5| image:: /dirac/images/math/5/d/1/5d1ccb648774ad117a597fb5cb7ad174.png
.. |image6| image:: /dirac/images/thumb/9/95/Dirac.png/41px-Dirac.png
.. |Powered by
MediaWiki| image:: /dirac/skins/common/images/poweredby_mediawiki_88x31.png
