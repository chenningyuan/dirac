Manual:ANALYZE:LOCALIZATION
===========================

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

Contents
--------

-  `1 \*LOCALIZATION -- Molecular orbital
   localization <#.2ALOCALIZATION_--_Molecular_orbital_localization>`_

   -  `1.1 .THRESH <#.THRESH>`_
   -  `1.2 .MAXITR <#.MAXITR>`_
   -  `1.3 .PRJLOC <#.PRJLOC>`_
   -  `1.4 .OWNBAS <#.OWNBAS>`_
   -  `1.5 .VECPRJ <#.VECPRJ>`_
   -  `1.6 .VECREF <#.VECREF>`_
   -  `1.7 .HESSUN <#.HESSUN>`_
   -  `1.8 .PRINT <#.PRINT>`_

\*LOCALIZATION -- Molecular orbital localization
================================================

Performs molecular orbitals localization using the Pipek-Mezey
criterion. The implementation uses an exponential parametrization and a
quasi-Newton minimization method instead of the Edminston-Ruedenberg
algorithm. For the moment, it only works with the *C*\ :sub:`1`\ 
symmetry and for the occupied orbitals.

.THRESH
-------

Convergence threshold.

*Default:*

::

    .THRESH
     1.0D-6

.MAXITR
-------

Maximum number of iterations.

*Default:*

::

    .MAXITR
     100

.PRJLOC
-------

Localization from projection analysis instead of Mulliken analysis.

::

    .PRJLOC
     1

*Default:*

::

    .PRJLOC
     0

.OWNBAS
-------

See `.OWNBAS </dirac/index.php/Manual:ANALYZE:PROJECTION#.OWNBAS>`_
under `\*PROJECTION </dirac/index.php/Manual:ANALYZE:PROJECTION>`_.

.VECPRJ
-------

See `.VECPRJ </dirac/index.php/Manual:ANALYZE:PROJECTION#.VECPRJ>`_
under `\*PROJECTION </dirac/index.php/Manual:ANALYZE:PROJECTION>`_.

.VECREF
-------

See `.VECREF </dirac/index.php/Manual:ANALYZE:PROJECTION#.VECREF>`_
under `\*PROJECTION </dirac/index.php/Manual:ANALYZE:PROJECTION>`_.

.HESSUN
-------

Hessian approximated by the unity matrix.

*Default:* Use the diagonal elements of the (diagonally dominant)
Hessian.

`|image0| </dirac/index.php/Image:HJAaJ.png>`_

**Programmers options**

.PRINT
------

Print level.

*Default:*

::

    .PRINT
     1

Retrieved from
"`http://wiki.chem.vu.nl/dirac/index.php/Manual:ANALYZE:LOCALIZATION <http://wiki.chem.vu.nl/dirac/index.php/Manual:ANALYZE:LOCALIZATION>`_"

Views
'''''

-  `Article </dirac/index.php/Manual:ANALYZE:LOCALIZATION>`_
-  `Discussion </dirac/index.php?title=Talk:Manual:ANALYZE:LOCALIZATION&action=edit>`_
-  `View
   source </dirac/index.php?title=Manual:ANALYZE:LOCALIZATION&action=edit>`_
-  `History </dirac/index.php?title=Manual:ANALYZE:LOCALIZATION&action=history>`_

Personal tools
''''''''''''''

-  `Log in / create
   account </dirac/index.php?title=Special:Userlogin&returnto=Manual:ANALYZE:LOCALIZATION>`_

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
   here </dirac/index.php/Special:Whatlinkshere/Manual:ANALYZE:LOCALIZATION>`_
-  `Related
   changes </dirac/index.php/Special:Recentchangeslinked/Manual:ANALYZE:LOCALIZATION>`_
-  `Upload file </dirac/index.php/Special:Upload>`_
-  `Special pages </dirac/index.php/Special:Specialpages>`_
-  `Printable
   version </dirac/index.php?title=Manual:ANALYZE:LOCALIZATION&printable=yes>`_
-  `Permanent
   link </dirac/index.php?title=Manual:ANALYZE:LOCALIZATION&oldid=3361>`_

`|Powered by MediaWiki| <http://www.mediawiki.org/>`_

-  This page was last modified 16:37, 9 March 2008.
-  This page has been accessed 1,060 times.
-  `Privacy policy </dirac/index.php/Diracwiki:Privacy_policy>`_
-  `About Diracwiki </dirac/index.php/Diracwiki:About>`_
-  `Disclaimers </dirac/index.php/Diracwiki:General_disclaimer>`_

.. |image0| image:: /dirac/images/thumb/6/65/HJAaJ.png/41px-HJAaJ.png
.. |Powered by
MediaWiki| image:: /dirac/skins/common/images/poweredby_mediawiki_88x31.png
