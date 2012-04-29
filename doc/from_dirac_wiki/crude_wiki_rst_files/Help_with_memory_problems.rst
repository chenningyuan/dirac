Help with memory problems
=========================

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

DIRAC will allocate the following memory segments:

-  Static allocation (94 MB; you can verify it using "make info")
-  WORK array for F77 "dynamic" allocation (specified by ./pam --mb)
-  Additional F90 dynamic allocations (depends on the calculation)

MEMGET errors
-------------

They look like this:

::

      --- SEVERE ERROR, PROGRAM WILL BE ABORTED ---

     Date and time (Linux) : Fri Nov 26 10:48:19 2010
     MEMGET ERROR, insufficient work space in memory
     Intentionally causing a crash with a Segmentation Fault

This means that you need to increase the WORK array (increase --mb).

MEMCHK errors
-------------

They look like this:

::

     MEMCHK ERROR, not a valid memget id in work(kfree-1)
     Text from calling routine : PSIDHF.DHFSCF (called from MEMREL)
     KFIRST,KFREE,IALLOC =         1         2         1
     found memory check :                        0
     expected           :               1234567890

This means that you have found a bug in the code. Please contact the
developers.

Other memory allocation errors
------------------------------

If you see a error message about memory problems in the queuing system
output (and not in the DIRAC output) it typically means that there is
not enough memory for F90 dynamic allocations.

You need to increase the memory limit that you specify to the queuing
system or perhaps decrease --mb (if possible).

Retrieved from
"`http://wiki.chem.vu.nl/dirac/index.php/Help\_with\_memory\_problems <http://wiki.chem.vu.nl/dirac/index.php/Help_with_memory_problems>`_"

Views
'''''

-  `Article </dirac/index.php/Help_with_memory_problems>`_
-  `Discussion </dirac/index.php?title=Talk:Help_with_memory_problems&action=edit>`_
-  `View
   source </dirac/index.php?title=Help_with_memory_problems&action=edit>`_
-  `History </dirac/index.php?title=Help_with_memory_problems&action=history>`_

Personal tools
''''''''''''''

-  `Log in / create
   account </dirac/index.php?title=Special:Userlogin&returnto=Help_with_memory_problems>`_

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
   here </dirac/index.php/Special:Whatlinkshere/Help_with_memory_problems>`_
-  `Related
   changes </dirac/index.php/Special:Recentchangeslinked/Help_with_memory_problems>`_
-  `Upload file </dirac/index.php/Special:Upload>`_
-  `Special pages </dirac/index.php/Special:Specialpages>`_
-  `Printable
   version </dirac/index.php?title=Help_with_memory_problems&printable=yes>`_
-  `Permanent
   link </dirac/index.php?title=Help_with_memory_problems&oldid=31481>`_

`|Powered by MediaWiki| <http://www.mediawiki.org/>`_

-  This page was last modified 15:02, 3 November 2011.
-  This page has been accessed 672 times.
-  `Privacy policy </dirac/index.php/Diracwiki:Privacy_policy>`_
-  `About Diracwiki </dirac/index.php/Diracwiki:About>`_
-  `Disclaimers </dirac/index.php/Diracwiki:General_disclaimer>`_

.. |Powered by
MediaWiki| image:: /dirac/skins/common/images/poweredby_mediawiki_88x31.png
