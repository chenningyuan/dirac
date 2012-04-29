Restarting
==========

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

A DIRAC calculation can be restarted at various stages depending on the
type of calculation. In general it should be enough to give the archive
file of a previous calculation as an argument to the ``pam`` script.
DIRAC will then find and use the available data for restart. In general
the program can restart from

-  ``DFCOEF`` Contains orbital coefficients, and will be used for
   restarting the SCF procedure.
-  ``AOPROPER`` Contains one electron integrals.
-  ``BSSMAT`` Contains four- to two-component transformation matrices,
   and allows you to skip the make\_bss step of a two component
   calculation (if you DIRAC version supports this).

It is also possible to restart a correlated calculation, see for example
the manual section of `MOLTRA </dirac/index.php/Manual:MOLTRA>`_, and
also the tips on `how to run a coupled-cluster calculation in separate
steps </dirac/index.php/How_to_run_a_calculation_in_separate_steps>`_.

Retrieved from
"`http://wiki.chem.vu.nl/dirac/index.php/Restarting <http://wiki.chem.vu.nl/dirac/index.php/Restarting>`_"

Views
'''''

-  `Article </dirac/index.php/Restarting>`_
-  `Discussion </dirac/index.php?title=Talk:Restarting&action=edit>`_
-  `View source </dirac/index.php?title=Restarting&action=edit>`_
-  `History </dirac/index.php?title=Restarting&action=history>`_

Personal tools
''''''''''''''

-  `Log in / create
   account </dirac/index.php?title=Special:Userlogin&returnto=Restarting>`_

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
   here </dirac/index.php/Special:Whatlinkshere/Restarting>`_
-  `Related
   changes </dirac/index.php/Special:Recentchangeslinked/Restarting>`_
-  `Upload file </dirac/index.php/Special:Upload>`_
-  `Special pages </dirac/index.php/Special:Specialpages>`_
-  `Printable
   version </dirac/index.php?title=Restarting&printable=yes>`_
-  `Permanent link </dirac/index.php?title=Restarting&oldid=16801>`_

`|Powered by MediaWiki| <http://www.mediawiki.org/>`_

-  This page was last modified 12:19, 21 April 2009.
-  This page has been accessed 975 times.
-  `Privacy policy </dirac/index.php/Diracwiki:Privacy_policy>`_
-  `About Diracwiki </dirac/index.php/Diracwiki:About>`_
-  `Disclaimers </dirac/index.php/Diracwiki:General_disclaimer>`_

.. |Powered by
MediaWiki| image:: /dirac/skins/common/images/poweredby_mediawiki_88x31.png
