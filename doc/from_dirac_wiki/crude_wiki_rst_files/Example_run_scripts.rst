Example run scripts
===================

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

Contents
--------

-  `1 Linux cluster "STALLO" at the University of Tromsø,
   Norway <#Linux_cluster_.22STALLO.22_at_the_University_of_Troms.C3.B8.2C_Norway>`_

Linux cluster "STALLO" at the University of Tromsø, Norway
----------------------------------------------------------

::

    #!/bin/bash

    # this is an example 8 processor job for 5 minutes
    # no commands before the PBS directives
    # please adapt [things in brackets]

    #PBS -S /bin/bash
    #PBS -A [your account in here]
    #PBS -l nodes=8:ppn=1:gige
    #PBS -l walltime=00:05:00
    #PBS -l pvmem=3800mb

    module load openmpi/1.3.2
    export LD_LIBRARY_PATH=/global/apps/intel/mkl/10.1.0.015/lib/em64t:$LD_LIBRARY_PATH
    export OMP_NUM_THREADS=1

    workdir=/global/work/$USER/$PBS_JOBID
    mkdir -p $workdir

    cd $PBS_O_WORKDIR
    /[edit path to pam]/pam -keepwrk -noarch -wrkdir $workdir -mpi 8 -mem 1700 -slavemem 1700 -aw 5 molecule.mol input.inp

    cd $PBS_O_WORKDIR
    rm -rf $workdir

    exit 0

Retrieved from
"`http://wiki.chem.vu.nl/dirac/index.php/Example\_run\_scripts <http://wiki.chem.vu.nl/dirac/index.php/Example_run_scripts>`_"

Views
'''''

-  `Article </dirac/index.php/Example_run_scripts>`_
-  `Discussion </dirac/index.php?title=Talk:Example_run_scripts&action=edit>`_
-  `View
   source </dirac/index.php?title=Example_run_scripts&action=edit>`_
-  `History </dirac/index.php?title=Example_run_scripts&action=history>`_

Personal tools
''''''''''''''

-  `Log in / create
   account </dirac/index.php?title=Special:Userlogin&returnto=Example_run_scripts>`_

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
   here </dirac/index.php/Special:Whatlinkshere/Example_run_scripts>`_
-  `Related
   changes </dirac/index.php/Special:Recentchangeslinked/Example_run_scripts>`_
-  `Upload file </dirac/index.php/Special:Upload>`_
-  `Special pages </dirac/index.php/Special:Specialpages>`_
-  `Printable
   version </dirac/index.php?title=Example_run_scripts&printable=yes>`_
-  `Permanent
   link </dirac/index.php?title=Example_run_scripts&oldid=17321>`_

`|Powered by MediaWiki| <http://www.mediawiki.org/>`_

-  This page was last modified 20:18, 9 June 2009.
-  This page has been accessed 1,755 times.
-  `Privacy policy </dirac/index.php/Diracwiki:Privacy_policy>`_
-  `About Diracwiki </dirac/index.php/Diracwiki:About>`_
-  `Disclaimers </dirac/index.php/Diracwiki:General_disclaimer>`_

.. |Powered by
MediaWiki| image:: /dirac/skins/common/images/poweredby_mediawiki_88x31.png
