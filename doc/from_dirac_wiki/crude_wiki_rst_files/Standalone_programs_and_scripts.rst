Standalone programs and scripts
===============================

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

DIRAC provides several standalone executable codes and scripts which are
either useful for us and perhaps also useful for you.

Contents
--------

-  `1 pamq.py script <#pamq.py_script>`_
-  `2 twofit.x executable <#twofit.x_executable>`_

pamq.py script
--------------

Script for submitting DIRAC for a specified queue system.

The current version supports few batch systems, like PBS, LoadLeveler,
and Grid Engine.

For example, save PBS as the default batch system, 1 node with 2 CPU's
(these are necessary parameters):

::

    pamq.py --config --qsys=pbs --nodes=1 --pnodes=2 

Submit and run DIRAC on 2 nodes in the 'express' queue (for 12 minutes).
One assumes that the queue system and number of CPU's per node are saved
in the configuration file:

::

    pamq.py --nodes=2 --time=00:12:00 --queue=express --pam_args="-noarch Mg.mol fsccd_IH.inp"

Only prepare script, do not sumbit

::

    pamq.py --nosub

For help, type

::

    pamq.py -h

twofit.x executable
-------------------

Serves for calculating spectroscopic parameters of diatomic molecules.

Retrieved from
"`http://wiki.chem.vu.nl/dirac/index.php/Standalone\_programs\_and\_scripts <http://wiki.chem.vu.nl/dirac/index.php/Standalone_programs_and_scripts>`_"

Views
'''''

-  `Article </dirac/index.php/Standalone_programs_and_scripts>`_
-  `Discussion </dirac/index.php?title=Talk:Standalone_programs_and_scripts&action=edit>`_
-  `View
   source </dirac/index.php?title=Standalone_programs_and_scripts&action=edit>`_
-  `History </dirac/index.php?title=Standalone_programs_and_scripts&action=history>`_

Personal tools
''''''''''''''

-  `Log in / create
   account </dirac/index.php?title=Special:Userlogin&returnto=Standalone_programs_and_scripts>`_

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
   here </dirac/index.php/Special:Whatlinkshere/Standalone_programs_and_scripts>`_
-  `Related
   changes </dirac/index.php/Special:Recentchangeslinked/Standalone_programs_and_scripts>`_
-  `Upload file </dirac/index.php/Special:Upload>`_
-  `Special pages </dirac/index.php/Special:Specialpages>`_
-  `Printable
   version </dirac/index.php?title=Standalone_programs_and_scripts&printable=yes>`_
-  `Permanent
   link </dirac/index.php?title=Standalone_programs_and_scripts&oldid=31571>`_

`|Powered by MediaWiki| <http://www.mediawiki.org/>`_

-  This page was last modified 15:13, 3 November 2011.
-  This page has been accessed 2,797 times.
-  `Privacy policy </dirac/index.php/Diracwiki:Privacy_policy>`_
-  `About Diracwiki </dirac/index.php/Diracwiki:About>`_
-  `Disclaimers </dirac/index.php/Diracwiki:General_disclaimer>`_

.. |Powered by
MediaWiki| image:: /dirac/skins/common/images/poweredby_mediawiki_88x31.png
