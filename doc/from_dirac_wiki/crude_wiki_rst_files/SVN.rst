SVN
===

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

Contents
--------

-  `1 Accessing the Dirac source code though
   SVN <#Accessing_the_Dirac_source_code_though_SVN>`_

   -  `1.1 About Subversion <#About_Subversion>`_

      -  `1.1.1 Main differences between CVS and
         SVN <#Main_differences_between_CVS_and_SVN>`_

   -  `1.2 Checking out the code <#Checking_out_the_code>`_
   -  `1.3 Updating the code <#Updating_the_code>`_
   -  `1.4 Adding, renaming, moving, and removing files and
      directories <#Adding.2C_renaming.2C_moving.2C_and_removing_files_and_directories>`_
   -  `1.5 Comparing code <#Comparing_code>`_
   -  `1.6 Commiting changes <#Commiting_changes>`_
   -  `1.7 Undoing changes <#Undoing_changes>`_
   -  `1.8 How to write useful commit
      messages <#How_to_write_useful_commit_messages>`_
   -  `1.9 Resolve conflicts <#Resolve_conflicts>`_
   -  `1.10 Good practices <#Good_practices>`_

Accessing the Dirac source code though SVN
==========================================

The DIRAC source is maintained using the Subversion (SVN) version
control system. An authorized developer (talk to HJJ) may use this tool
to check out and make changes to DIRAC.

About Subversion
----------------

`Subversion <http://subversion.tigris.org/>`_ is a program that helps
multiple developers work together on a program simultaneously. Every
revision to the source code is retained, so if something goes wrong it
is always possible to go back to a previous state. SVN comes with an
excellent online (and print)
`manual <http://svnbook.red-bean.com/en/1.0/>`_.

Use also the command-line help

::

    svn help

Main differences between CVS and SVN
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SVN was created in an attempt to fix many of the fundamental problems of
CVS. It mostly succeeded - it is very difficult to find someone who
liked CVS better. See the SVN manual chapter on `SVN for CVS
users <http://svnbook.red-bean.com/en/1.0/apa.html>`_.

-  You can move files and directories (``svn mv``) also across
   directories by keeping the history.
-  You can remove files and directories (``svn rm``) of course by
   keeping the history.
-  Branches and tags are handled slightly differently, see the
   `manual <http://svnbook.red-bean.com/en/1.0/ch04.html>`_ for more
   information.
-  Revision numbers are not per-file but label the state of the
   repository, commits are atomic (either entire commit is ok or no
   commit takes place).
-  SVN makes it more difficult to commit files in conflict.
-  Don't have to do ``svn up`` if you rather wanted to see the status of
   changes (``svn st``).

Checking out the code
---------------------

To obtain a copy of the latest DIRAC source tree you can do

::

    svn co http://username@dirac.blah.blah/path/to/svn/trunk

The word ``co`` is short for "check out". The final directory ``trunk``
indicates that we want a copy of the main development version.

When the "username" for the remote dirac repository is different from
your login name on the connecting client, please do the check out as

::

    svn --username <your_dirac_repository_name> co http://dirac.blah.blah/path/to/svn/trunk <yourDirac>

where <yourDirac> is the name of your local DIRAC directory. Leaving out
<yourDirac> in the command the directory will have same name as in the
repository - ``trunk``. You may choose the name <yourDirac> according to
characteristics of your local version, for example
"dirac\_gfortran\_gdb", "dirac\_ifort\_optim\_serial" etc...

 For more see the `manual <http://svnbook.red-bean.com/en/1.0/>`_ for
information about tags and branches.

Updating the code
-----------------

To update the entire directory (with its subdirectories) use

::

    svn up

or

::

    svn update

You can update only a specific file or subdirectory and also update to a
specific revision number or date

::

    svn up some_file.F90 --revision {2008-04-25}

To update the whole DIRAC program suite you must go the the main
directory.

Adding, renaming, moving, and removing files and directories
------------------------------------------------------------

That's right - with SVN you can not only add files and directories

::

    svn add some_directory
    svn add file1 file2.ps file3.F90

but also move (or rename) files and directories

::

    svn mv file1 better_name_for_file1
    svn mv directory1/some_file directory2/better_name

by keeping the history.

You can remove files and directories

::

    svn rm obsolete_directory

and you don't need to delete them first before removing them.

After you have added, moved, or deleted files, you need to commit the
changes with

::

    svn ci

(see also further below).

Comparing code
--------------

To compare a file against its most recent *local* version use

::

    svn diff some_file_name.F90

A plain

::

    svn diff

will diff the entire directory.

You can also compare with a specific revision number

::

    svn diff some_file.F90 -r20

or date

::

    svn diff some_file.F90 --revision {2008-04-25}

Compare two different revisions (N and M) of the same file

::

    svn diff -r N:M some_file.F90 

You can list all files that you have changed after your last update with

::

    svn st

or

::

    svn status

A quite useful command to find out who changed a specific line of code
when exactly is

::

    svn blame some_file

(not to be used to blame developers).

Commiting changes
-----------------

Before checking in new code, update first

::

    svn up

then commit

::

    svn ci some_file

or

::

    svn ci

to commit all modifications that you see with `` svn st `` and write a
useful commit message (see below) with your favorite editor (set by the
environment variable ``$EDITOR``).

You can also tell SVN to get the message from a file by passing the
filename with the --file switch:

::

    svn ci --file logmsg

If you're in your editor writing a commit message and decide that you
want to cancel your commit, you can just quit your editor without saving
changes. If you've already saved your commit message, simply delete the
text and save again.

Undoing changes
---------------

Suppose you discovered that the changeset (long time ago) between
revisions 415 and 416 was a mistake. Just to check we do a diff

::

    svn diff -r 415:416 file

and this we can undo with the command

::

    svn merge -r 416:415 file

It is possible that you get conflicts which you have to resolve (see
further below).

How to write useful commit messages
-----------------------------------

First look at old commit messages with

::

    svn log

Useful commit messages serve as code development documentation, make it
easy for people to search in the future.

Useful commit messages keep the developers informed about bugfixes and
architectural changes. Since developers often use the latest revision
for calculations, it is good to be more specific than just writing
"bugfix" so other developers know whether they have to worry about their
results. Other developers can learn from useful bugfix commit logs.

It is not necessary to write your name/initials or list the modified
files, SVN does this for you.

Resolve conflicts
-----------------

You may observe a conflict after an update (``svn up``)

::

    C    your_modified_file.F90
    Updated to revision 23.

This means that the changes from the server overlapped with your own,
and now you have to manually choose between them.

You need to do one of three things:

-  Merge the conflicted text “by hand” (by examining and editing the
   conflict markers within the file).
-  Copy one of the temporary files (created by SVN) on top of your
   working file.
-  Run ``svn revert your_modified_file.F90`` to throw away all of your
   local changes.

Once you've resolved the conflict, you need to let SVN know by running

::

    svn resolved your_modified_file.F90

Good practices
--------------

-  Don't break the trunk.

-  If it builds on your machine the changes still can break the trunk.
   Two favorite mistakes are forgetting to add files and misplacing CPP
   flags (it can build with VAR\_MPI set but not without and vice
   versa). It is a good practice to keep separate working copies of the
   repository where you don't edit, only update - ideally one for a
   sequential and one for a parallel build. After your commit, try to
   update the testing working copies and check whether it builds
   correctly there.

-  Commit changes in coherent units (don't commit unrelated changes in
   single commits). This is because if you commit something that is a
   good idea and something which is not in one single step, it is
   cumbersome to rollback the bad idea without undoing also the good
   one. It is also easier to collect patches for releases and updates
   based on coherent commits.

Retrieved from
"`http://wiki.chem.vu.nl/dirac/index.php/SVN <http://wiki.chem.vu.nl/dirac/index.php/SVN>`_"

Views
'''''

-  `Article </dirac/index.php/SVN>`_
-  `Discussion </dirac/index.php/Talk:SVN>`_
-  `View source </dirac/index.php?title=SVN&action=edit>`_
-  `History </dirac/index.php?title=SVN&action=history>`_

Personal tools
''''''''''''''

-  `Log in / create
   account </dirac/index.php?title=Special:Userlogin&returnto=SVN>`_

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

-  `What links here </dirac/index.php/Special:Whatlinkshere/SVN>`_
-  `Related changes </dirac/index.php/Special:Recentchangeslinked/SVN>`_
-  `Upload file </dirac/index.php/Special:Upload>`_
-  `Special pages </dirac/index.php/Special:Specialpages>`_
-  `Printable version </dirac/index.php?title=SVN&printable=yes>`_
-  `Permanent link </dirac/index.php?title=SVN&oldid=14361>`_

`|Powered by MediaWiki| <http://www.mediawiki.org/>`_

-  This page was last modified 12:44, 3 November 2008.
-  This page has been accessed 1,772 times.
-  `Privacy policy </dirac/index.php/Diracwiki:Privacy_policy>`_
-  `About Diracwiki </dirac/index.php/Diracwiki:About>`_
-  `Disclaimers </dirac/index.php/Diracwiki:General_disclaimer>`_

.. |Powered by
MediaWiki| image:: /dirac/skins/common/images/poweredby_mediawiki_88x31.png
