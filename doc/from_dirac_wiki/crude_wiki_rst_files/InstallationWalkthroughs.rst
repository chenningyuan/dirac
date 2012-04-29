InstallationWalkthroughs
========================

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

Hyperion (2012-02-13)
^^^^^^^^^^^^^^^^^^^^^

First set the correct environment. For this unload all possibly
conflicting modules:

::

    module unload default-intel-cc default-intel-fc default-intel-mkl default-mpt-sgi

Then load:

::

    module load intel-cc-12/12.0.4 intel-fc-12/12.0.4 intel-mkl/12.0.4

Also hook up to own 64bit integer OpenMPI build (ask Radovan for correct
path):

::

    export PATH=$PATH':/path/to/openmpi/bin'
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH':/path/to/openmpi/lib'

::

    export MATH_ROOT='/divers/intel/Compiler/12.0.4/mkl'
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH':/divers/intel/Compiler/12.0.4/mkl/lib/intel64'

Now setup DIRAC:

::

    ./setup --fc=mpif90 --cc=mpicc --cxx=icpc --int64 -D ENABLE_LARGE_MOLECULES=1

Finally build the code:

::

    cd build
    make -j16

Retrieved from
"`http://wiki.chem.vu.nl/dirac/index.php/InstallationWalkthroughs <http://wiki.chem.vu.nl/dirac/index.php/InstallationWalkthroughs>`_"

Views
'''''

-  `Article </dirac/index.php/InstallationWalkthroughs>`_
-  `Discussion </dirac/index.php?title=Talk:InstallationWalkthroughs&action=edit>`_
-  `View
   source </dirac/index.php?title=InstallationWalkthroughs&action=edit>`_
-  `History </dirac/index.php?title=InstallationWalkthroughs&action=history>`_

Personal tools
''''''''''''''

-  `Log in / create
   account </dirac/index.php?title=Special:Userlogin&returnto=InstallationWalkthroughs>`_

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
   here </dirac/index.php/Special:Whatlinkshere/InstallationWalkthroughs>`_
-  `Related
   changes </dirac/index.php/Special:Recentchangeslinked/InstallationWalkthroughs>`_
-  `Upload file </dirac/index.php/Special:Upload>`_
-  `Special pages </dirac/index.php/Special:Specialpages>`_
-  `Printable
   version </dirac/index.php?title=InstallationWalkthroughs&printable=yes>`_
-  `Permanent
   link </dirac/index.php?title=InstallationWalkthroughs&oldid=33741>`_

`|Powered by MediaWiki| <http://www.mediawiki.org/>`_

-  This page was last modified 11:22, 13 February 2012.
-  This page has been accessed 250 times.
-  `Privacy policy </dirac/index.php/Diracwiki:Privacy_policy>`_
-  `About Diracwiki </dirac/index.php/Diracwiki:About>`_
-  `Disclaimers </dirac/index.php/Diracwiki:General_disclaimer>`_

.. |Powered by
MediaWiki| image:: /dirac/skins/common/images/poweredby_mediawiki_88x31.png
