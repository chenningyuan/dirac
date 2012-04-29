How to safely move your uncommited files to SVN
===============================================

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

If you decide/need to move your modified files to SVN **without**
committing them to CVS first, this is a safe way:

-  Update your CVS trunk working copy or branch working copy
-  Check out a fresh SVN working copy (either trunk or a specific
   branch)
-  Move modified files by hand to the new SVN working copy
-  **Important**: Check with ``svn status`` whether SVN marks the same
   files as modified (the list of files given by ``svn status`` should
   look like the list of files given by ``cvs update``). If this is not
   the case and/or you have doubts, don't delete the CVS working copy
   too early.

Retrieved from
"`http://wiki.chem.vu.nl/dirac/index.php/How\_to\_safely\_move\_your\_uncommited\_files\_to\_SVN <http://wiki.chem.vu.nl/dirac/index.php/How_to_safely_move_your_uncommited_files_to_SVN>`_"

Views
'''''

-  `Article </dirac/index.php/How_to_safely_move_your_uncommited_files_to_SVN>`_
-  `Discussion </dirac/index.php?title=Talk:How_to_safely_move_your_uncommited_files_to_SVN&action=edit>`_
-  `View
   source </dirac/index.php?title=How_to_safely_move_your_uncommited_files_to_SVN&action=edit>`_
-  `History </dirac/index.php?title=How_to_safely_move_your_uncommited_files_to_SVN&action=history>`_

Personal tools
''''''''''''''

-  `Log in / create
   account </dirac/index.php?title=Special:Userlogin&returnto=How_to_safely_move_your_uncommited_files_to_SVN>`_

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
   here </dirac/index.php/Special:Whatlinkshere/How_to_safely_move_your_uncommited_files_to_SVN>`_
-  `Related
   changes </dirac/index.php/Special:Recentchangeslinked/How_to_safely_move_your_uncommited_files_to_SVN>`_
-  `Upload file </dirac/index.php/Special:Upload>`_
-  `Special pages </dirac/index.php/Special:Specialpages>`_
-  `Printable
   version </dirac/index.php?title=How_to_safely_move_your_uncommited_files_to_SVN&printable=yes>`_
-  `Permanent
   link </dirac/index.php?title=How_to_safely_move_your_uncommited_files_to_SVN&oldid=6731>`_

`|Powered by MediaWiki| <http://www.mediawiki.org/>`_

-  This page was last modified 09:57, 16 May 2008.
-  This page has been accessed 250 times.
-  `Privacy policy </dirac/index.php/Diracwiki:Privacy_policy>`_
-  `About Diracwiki </dirac/index.php/Diracwiki:About>`_
-  `Disclaimers </dirac/index.php/Diracwiki:General_disclaimer>`_

.. |Powered by
MediaWiki| image:: /dirac/skins/common/images/poweredby_mediawiki_88x31.png
