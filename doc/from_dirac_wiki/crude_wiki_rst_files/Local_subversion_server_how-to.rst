Local subversion server how-to
==============================

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

This script makes only sense when interfacing a local SVN server with
the central CVS server. Synchronization of a local SVN server with
another SVN server is easier and I will replace this how-to as soon as
we have moved to SVN.

 This script will set up a new subversion repository for DIRAC together
with a sync working copy to sync the repository with the cvs Odense.

It is useful to set:

::

    alias cvs='cvs -I .svn'

::

    #!/bin/bash

    # edit from here ...

      # WARNING: make sure the following directory exists:
      dirac_directory=/home/bast/dirac

      # WARNING: the following directory will be deleted if it exists already:
      dirac_sync_directory=sync

      # WARNING: the following repository will be wiped out if it exists already:
      subversion_repository=/home/bast/svn_dirac

      will_be_deleted=/home/bast/delete_me

      # contains files that svn will not version (ignore):
      what_to_ignore=$PWD/what_to_ignore

      # adapt it to your needs:
    cat << EOF > $what_to_ignore
    .cvsignore
    CVS
    *.o
    *.i
    Makefile*
    include
    EOF

    # ... to here

    # cvs checkout to sync directory

    cd $dirac_directory
    rm -rf $dirac_sync_directory
    cvs co -d $dirac_sync_directory Dirac

    # initial svn repository layout

    umask=077

    rm -rf $subversion_repository
    mkdir $subversion_repository

    svnadmin create --fs-type fsfs $subversion_repository

    rm -rf $will_be_deleted
    mkdir $will_be_deleted
    cd $will_be_deleted

    mkdir dirac
    mkdir dirac/trunk
    mkdir dirac/branches
    mkdir dirac/tags

    svn import . file://$subversion_repository -m 'initial repository layout'

    rm -rf $will_be_deleted/*
    cd $will_be_deleted

    svn co file://$subversion_repository/dirac/trunk dirac

    # copy dirac files to $will_be_deleted and remove cvs related files for initial commit to svn

    cp -r $dirac_directory/$dirac_sync_directory/* $will_be_deleted/dirac
    cd $will_be_deleted/dirac
    find . -name ".cvsignore" -exec rm '{}' \;
    find . -name "CVS" -exec rm -rf '{}' &> /dev/null \;

    # initial commit to svn

    for file in *; do
        svn add $file
    done
    svn commit -m "initial commit"

    # set ignores

    svn update
    find . -type d | grep -v '.svn' | grep -v CVS | xargs svn propset svn:ignore -F $what_to_ignore
    svn commit -m "ignores set"

    # move svn related files to sync directory

    tar -cf $dirac_directory/$dirac_sync_directory/untar_me.tar `find -type d | grep ".svn"`
    cd $dirac_directory/$dirac_sync_directory
    tar -xf untar_me.tar >& /dev/null

    # cleanup

    rm -rf $will_be_deleted
    rm -f untar_me.tar
    rm -f $what_to_ignore

    exit 0

After this is done you can update your svn repository with recent
changes in cvs Odense:

::

    cvs update; svn commit

... or sync the cvs Odense with your changes:

::

    svn update; cvs commit

Retrieved from
"`http://wiki.chem.vu.nl/dirac/index.php/Local\_subversion\_server\_how-to <http://wiki.chem.vu.nl/dirac/index.php/Local_subversion_server_how-to>`_"

Views
'''''

-  `Article </dirac/index.php/Local_subversion_server_how-to>`_
-  `Discussion </dirac/index.php?title=Talk:Local_subversion_server_how-to&action=edit>`_
-  `View
   source </dirac/index.php?title=Local_subversion_server_how-to&action=edit>`_
-  `History </dirac/index.php?title=Local_subversion_server_how-to&action=history>`_

Personal tools
''''''''''''''

-  `Log in / create
   account </dirac/index.php?title=Special:Userlogin&returnto=Local_subversion_server_how-to>`_

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
   here </dirac/index.php/Special:Whatlinkshere/Local_subversion_server_how-to>`_
-  `Related
   changes </dirac/index.php/Special:Recentchangeslinked/Local_subversion_server_how-to>`_
-  `Upload file </dirac/index.php/Special:Upload>`_
-  `Special pages </dirac/index.php/Special:Specialpages>`_
-  `Printable
   version </dirac/index.php?title=Local_subversion_server_how-to&printable=yes>`_
-  `Permanent
   link </dirac/index.php?title=Local_subversion_server_how-to&oldid=6061>`_

`|Powered by MediaWiki| <http://www.mediawiki.org/>`_

-  This page was last modified 08:05, 14 May 2008.
-  This page has been accessed 115 times.
-  `Privacy policy </dirac/index.php/Diracwiki:Privacy_policy>`_
-  `About Diracwiki </dirac/index.php/Diracwiki:About>`_
-  `Disclaimers </dirac/index.php/Diracwiki:General_disclaimer>`_

.. |Powered by
MediaWiki| image:: /dirac/skins/common/images/poweredby_mediawiki_88x31.png
