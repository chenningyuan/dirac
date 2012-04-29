Manual:WAVE FUNCTION:RESOLVE
============================

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

Contents
--------

-  `1 \*RESOLVE -- Resolve open-shell
   states <#.2ARESOLVE_--_Resolve_open-shell_states>`_

   -  `1.1 .PRINT <#.PRINT>`_
   -  `1.2 .INTFLG <#.INTFLG>`_
   -  `1.3 .SCREEN <#.SCREEN>`_
   -  `1.4 .SCHEME <#.SCHEME>`_

\*RESOLVE -- Resolve open-shell states
======================================

The individual electronic states of an average-of-configurations
calculation are resolved by a full CI in the small open-shell
manifold(s) using the `GOSCIP module </dirac/index.php/Manual:GOSCIP>`_.

Note that even if this is a small CI diagonalization one needs in
principle to generate *all* AO integrals for the initial 4-index
transformation. The cost of calculation can be greatly reduced invoking
integral screening (controlled by `.SCREEN <#.SCREEN>`_) in the 4-index
transformation, possibly also by leaving out at least SS integrals
(using `.INTFLG <#.INTFLG>`_).

.PRINT
------

Print level.

*Default:*

::

    .PRINT
     0

.INTFLG
-------

Specify what two-electron integrals to include (see
`.INTFLG </dirac/index.php/Manual:HAMILTONIAN#.INTFLG>`_ under
`\*\*HAMILTONIAN </dirac/index.php/Manual:HAMILTONIAN>`_).

*Default:* `.INTFLG </dirac/index.php/Manual:HAMILTONIAN#.INTFLG>`_ from
`\*\*HAMILTONIAN </dirac/index.php/Manual:HAMILTONIAN>`_.

.SCREEN
-------

Screening threshold in the 4-index transformation.

A negative number deactivates screening.

*Default:*

::

    .SCREEN
     1.0D-14

.SCHEME
-------

Integral transformation algorithm.\ :sup:`[*please
`clarify </dirac/index.php/Clarify>`_*]`\ 

*Default:*

::

    .SCHEME
     4

Retrieved from
"`http://wiki.chem.vu.nl/dirac/index.php/Manual:WAVE\_FUNCTION:RESOLVE <http://wiki.chem.vu.nl/dirac/index.php/Manual:WAVE_FUNCTION:RESOLVE>`_"

`Category </dirac/index.php/Special:Categories>`_: `Articles needing
clarification </dirac/index.php/Category:Articles_needing_clarification>`_

Views
'''''

-  `Article </dirac/index.php/Manual:WAVE_FUNCTION:RESOLVE>`_
-  `Discussion </dirac/index.php?title=Talk:Manual:WAVE_FUNCTION:RESOLVE&action=edit>`_
-  `View
   source </dirac/index.php?title=Manual:WAVE_FUNCTION:RESOLVE&action=edit>`_
-  `History </dirac/index.php?title=Manual:WAVE_FUNCTION:RESOLVE&action=history>`_

Personal tools
''''''''''''''

-  `Log in / create
   account </dirac/index.php?title=Special:Userlogin&returnto=Manual:WAVE_FUNCTION:RESOLVE>`_

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
   here </dirac/index.php/Special:Whatlinkshere/Manual:WAVE_FUNCTION:RESOLVE>`_
-  `Related
   changes </dirac/index.php/Special:Recentchangeslinked/Manual:WAVE_FUNCTION:RESOLVE>`_
-  `Upload file </dirac/index.php/Special:Upload>`_
-  `Special pages </dirac/index.php/Special:Specialpages>`_
-  `Printable
   version </dirac/index.php?title=Manual:WAVE_FUNCTION:RESOLVE&printable=yes>`_
-  `Permanent
   link </dirac/index.php?title=Manual:WAVE_FUNCTION:RESOLVE&oldid=11371>`_

`|Powered by MediaWiki| <http://www.mediawiki.org/>`_

-  This page was last modified 11:03, 16 July 2008.
-  This page has been accessed 1,082 times.
-  `Privacy policy </dirac/index.php/Diracwiki:Privacy_policy>`_
-  `About Diracwiki </dirac/index.php/Diracwiki:About>`_
-  `Disclaimers </dirac/index.php/Diracwiki:General_disclaimer>`_

.. |Powered by
MediaWiki| image:: /dirac/skins/common/images/poweredby_mediawiki_88x31.png
