Manual:ANALYZE:RHO1
===================

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

Contents
--------

-  `1 \*RHO1 -- Print density along lines between atomic
   centers <#.2ARHO1_--_Print_density_along_lines_between_atomic_centers>`_

   -  `1.1 .PRINT <#.PRINT>`_
   -  `1.2 .MESH <#.MESH>`_

\*RHO1 -- Print density along lines between atomic centers
==========================================================

For each unique pair (A,B) of centers a formatted file is created,
containing four columns of numbers: The origin is placed at atom A, and
the first column gives the coordinate in Angstroms along the line
between centers A and B. The following columns contain the total, large
and small component density, respectively (in atomic units). The file is
named 'RH',NAMN(ICENTA),NAMN(ICENTB),IDEGB where ICENTB is the name of
symmetry-independent center B and IDEGB the numbering of the dependent
center. Underscores are inserted wherever there are blanks. The
generation of the density requires the coefficient file DFCOEF. Note
that one can in principle obtain the density along any line simply by
introducing ghost centers.

For uranyl one may use the command

::

    pam -incmo -get ``RHU___O____1  RHU___O____2 ...

.PRINT
------

Print level.

*Default:*

::

    .PRINT
     0

.MESH
-----

Step length for grid.

*Default:*

::

    .MESH
     0.01

Retrieved from
"`http://wiki.chem.vu.nl/dirac/index.php/Manual:ANALYZE:RHO1 <http://wiki.chem.vu.nl/dirac/index.php/Manual:ANALYZE:RHO1>`_"

Views
'''''

-  `Article </dirac/index.php/Manual:ANALYZE:RHO1>`_
-  `Discussion </dirac/index.php?title=Talk:Manual:ANALYZE:RHO1&action=edit>`_
-  `View
   source </dirac/index.php?title=Manual:ANALYZE:RHO1&action=edit>`_
-  `History </dirac/index.php?title=Manual:ANALYZE:RHO1&action=history>`_

Personal tools
''''''''''''''

-  `Log in / create
   account </dirac/index.php?title=Special:Userlogin&returnto=Manual:ANALYZE:RHO1>`_

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
   here </dirac/index.php/Special:Whatlinkshere/Manual:ANALYZE:RHO1>`_
-  `Related
   changes </dirac/index.php/Special:Recentchangeslinked/Manual:ANALYZE:RHO1>`_
-  `Upload file </dirac/index.php/Special:Upload>`_
-  `Special pages </dirac/index.php/Special:Specialpages>`_
-  `Printable
   version </dirac/index.php?title=Manual:ANALYZE:RHO1&printable=yes>`_
-  `Permanent
   link </dirac/index.php?title=Manual:ANALYZE:RHO1&oldid=3711>`_

`|Powered by MediaWiki| <http://www.mediawiki.org/>`_

-  This page was last modified 09:45, 10 March 2008.
-  This page has been accessed 736 times.
-  `Privacy policy </dirac/index.php/Diracwiki:Privacy_policy>`_
-  `About Diracwiki </dirac/index.php/Diracwiki:About>`_
-  `Disclaimers </dirac/index.php/Diracwiki:General_disclaimer>`_

.. |Powered by
MediaWiki| image:: /dirac/skins/common/images/poweredby_mediawiki_88x31.png
