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

