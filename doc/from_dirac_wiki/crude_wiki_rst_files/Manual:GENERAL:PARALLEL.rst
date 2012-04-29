Manual:GENERAL:PARALLEL
=======================

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

Contents
--------

-  `1 \*PARALLEL -- Parallelization
   directives <#.2APARALLEL_--_Parallelization_directives>`_

   -  `1.1 .PRINT <#.PRINT>`_
   -  `1.2 .NTASK <#.NTASK>`_

\*PARALLEL -- Parallelization directives
========================================

`|image0| </dirac/index.php/Image:HJAaJ.png>`_

**Programmers options**

.PRINT
------

Print level for the parallel calculation.

*Default:*

::

    .PRINT
     0

A print level of at least 2 is needed in order to be able to evaluate
the parallelization efficiency. A complete timing for all nodes will be
given if the print level is 4 or higher.

.NTASK
------

Number of tasks to send to each node when distributing the calculation
of two-electron integrals.

*Default:*

::

    .NTASK
     1

A task is defined as a shell of atomic integrals, a shell being an input
block. One may therefore increase the number of shells given to each
node in order to reduce the amount of communication. However, the
program uses dynamical allocation of work to each node, and thus this
option should be used with some care, as too large tasks may cause the
dynamical load balancing to fail, giving an overall decrease in
efficiency. The parallelization is also very coarse grained, so that the
amount of communication seldom represents any significant problem.

Retrieved from
"`http://wiki.chem.vu.nl/dirac/index.php/Manual:GENERAL:PARALLEL <http://wiki.chem.vu.nl/dirac/index.php/Manual:GENERAL:PARALLEL>`_"

Views
'''''

-  `Article </dirac/index.php/Manual:GENERAL:PARALLEL>`_
-  `Discussion </dirac/index.php?title=Talk:Manual:GENERAL:PARALLEL&action=edit>`_
-  `View
   source </dirac/index.php?title=Manual:GENERAL:PARALLEL&action=edit>`_
-  `History </dirac/index.php?title=Manual:GENERAL:PARALLEL&action=history>`_

Personal tools
''''''''''''''

-  `Log in / create
   account </dirac/index.php?title=Special:Userlogin&returnto=Manual:GENERAL:PARALLEL>`_

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
   here </dirac/index.php/Special:Whatlinkshere/Manual:GENERAL:PARALLEL>`_
-  `Related
   changes </dirac/index.php/Special:Recentchangeslinked/Manual:GENERAL:PARALLEL>`_
-  `Upload file </dirac/index.php/Special:Upload>`_
-  `Special pages </dirac/index.php/Special:Specialpages>`_
-  `Printable
   version </dirac/index.php?title=Manual:GENERAL:PARALLEL&printable=yes>`_
-  `Permanent
   link </dirac/index.php?title=Manual:GENERAL:PARALLEL&oldid=3231>`_

`|Powered by MediaWiki| <http://www.mediawiki.org/>`_

-  This page was last modified 16:32, 9 March 2008.
-  This page has been accessed 1,010 times.
-  `Privacy policy </dirac/index.php/Diracwiki:Privacy_policy>`_
-  `About Diracwiki </dirac/index.php/Diracwiki:About>`_
-  `Disclaimers </dirac/index.php/Diracwiki:General_disclaimer>`_

.. |image0| image:: /dirac/images/thumb/6/65/HJAaJ.png/41px-HJAaJ.png
.. |Powered by
MediaWiki| image:: /dirac/skins/common/images/poweredby_mediawiki_88x31.png
