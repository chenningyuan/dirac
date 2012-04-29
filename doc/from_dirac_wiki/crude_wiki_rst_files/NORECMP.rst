NORECMP
=======

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

Contents
--------

-  `1 &CCSORT -- Sorting module <#.26CCSORT_--_Sorting_module>`_

   -  `1.1 NORECMP <#NORECMP>`_

&CCSORT -- Sorting module
=========================

NORECMP
-------

For the open shell calculation, turn off this flag. Use it with USEOE
option.

::

    &RELCCSD NELEC=4,3 TIMING=T, IPRNT=1, DOSORT=T, DOENER=T &END
    &CCSORT USEOE=F, NORECMP=F &END
    &CCENER DOMP2=T, DOCCSD=T, DOCCSDT=F &END

Retrieved from
"`http://wiki.chem.vu.nl/dirac/index.php/NORECMP <http://wiki.chem.vu.nl/dirac/index.php/NORECMP>`_"

Views
'''''

-  `Article </dirac/index.php/NORECMP>`_
-  `Discussion </dirac/index.php?title=Talk:NORECMP&action=edit>`_
-  `View source </dirac/index.php?title=NORECMP&action=edit>`_
-  `History </dirac/index.php?title=NORECMP&action=history>`_

Personal tools
''''''''''''''

-  `Log in / create
   account </dirac/index.php?title=Special:Userlogin&returnto=NORECMP>`_

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

-  `What links here </dirac/index.php/Special:Whatlinkshere/NORECMP>`_
-  `Related
   changes </dirac/index.php/Special:Recentchangeslinked/NORECMP>`_
-  `Upload file </dirac/index.php/Special:Upload>`_
-  `Special pages </dirac/index.php/Special:Specialpages>`_
-  `Printable version </dirac/index.php?title=NORECMP&printable=yes>`_
-  `Permanent link </dirac/index.php?title=NORECMP&oldid=20521>`_

`|Powered by MediaWiki| <http://www.mediawiki.org/>`_

-  This page was last modified 12:52, 19 May 2010.
-  This page has been accessed 100 times.
-  `Privacy policy </dirac/index.php/Diracwiki:Privacy_policy>`_
-  `About Diracwiki </dirac/index.php/Diracwiki:About>`_
-  `Disclaimers </dirac/index.php/Diracwiki:General_disclaimer>`_

.. |Powered by
MediaWiki| image:: /dirac/skins/common/images/poweredby_mediawiki_88x31.png
