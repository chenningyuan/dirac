CVS2SVN
=======

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

Contents
--------

-  `1 Basic conversion <#Basic_conversion>`_
-  `2 Conversion via svndumpfilter <#Conversion_via_svndumpfilter>`_
-  `3 Ignores <#Ignores>`_
-  `4 Backup <#Backup>`_
-  `5 Commit hooks <#Commit_hooks>`_

Basic conversion
----------------

Set the variables $SVN\_DIRECTORY and $CVS\_DIRECTORY.

First remove unwanted directories directly in the $CVS\_DIRECTORY. These
are CVSROOT and www/passwd. It is difficult to remove them later during
the conversion.

Then

::

    cvs2svn -s $SVN_DIRECTORY $CVS_DIRECTORY

Will create a new SVN repository which contains all old branches and
tags. The process can take 30 to 80 minutes for nearly 7000 revisions.

Look at the results with a command like

::

    svn ls file://$SVN_DIRECTORY

I get

-  branches/

   -  Dirac04/
   -  Dirac07/
   -  LV\_Gaunt/
      ..

-  tags/
   ..
-  trunk/

To check out the trunk do

::

    svn co file://$SVN_DIRECTORY/trunk

Conversion via svndumpfilter
----------------------------

**Problem**: This route described below doesn't work for me, during the
svnadmin load I get errors of the type:

::

    svnadmin: File not found: transaction '771-1', path 'branches/unlabeled-1.8.6/CVSROOT'
        * deleting path : branches/unlabeled-1.8.6/CVSROOT ...

No idea what the unlabeled branches are. I think we need to remove the
stuff (trunk/CVSROOT and trunk/www/passwd) the hard way from the CVS
repository before the conversion.

The idea was to have a conversion which would filter out unwanted
directories (trunk/CVSROOT and trunk/www/passwd) before the SVN
repository is created.

Set the variables $SVN\_DIRECTORY and $CVS\_DIRECTORY, $DUMP\_FILE and
$FILTERED\_DUMP\_FILE.

First we convert the CVS repository to a SVN dump file:

::

    cvs2svn --dumpfile=$DUMP_FILE $CVS_DIRECTORY

Then we filter out the unwanted directories:

::

    cat $DUMP_FILE | svndumpfilter exclude trunk/CVSROOT    --drop-empty-revs --renumber-revs \
                   | svndumpfilter exclude trunk/www/passwd --drop-empty-revs --renumber-revs > $FILTERED_DUMP_FILE

Next we create an empty SVN repository:

::

    umask=077
    svnadmin create --fs-type fsfs $SVN_DIRECTORY

And finally load the filtered dump file:

::

    svnadmin load $SVN_DIRECTORY < $FILTERED_DUMP_FILE

Ignores
-------

Not sure we really need it - but this is a way to add the ignore
property to directories and subdirectories.

Create a file which contains ignore patterns (set $IGNORE\_FILE):

::

    cat << EOF > $IGNORE_FILE
    Makefile
    Makefile.config
    *.o
    *.i
    *.old
    *.depend
    *.mod
    EOF

Then add the ignore property to the current directory (".") and in this
case to all subdirectories:

::

    find . -type d | grep -v '.svn' | xargs svn propset svn:ignore -F $IGNORE_FILE

For a specific directory only, use:

::

    svn propset svn:ignore -F $IGNORE_FILE /this/directory/only

Finally commit the change:

::

    svn ci

Backup
------

An easy way to backup the entire repository into one single file is:

::

    svnadmin dump $SVN_DIRECTORY > $BACKUP_FILE

This is conveniently done with as a crontab job.

The repository can be restored with the following two commands:

::

    svnadmin create --fs-type fsfs $SVN_DIRECTORY
    svnadmin load $SVN_DIRECTORY < $BACKUP_FILE

Commit hooks
------------

Download
`commit-email.pl <http://svn.collab.net/repos/svn/trunk/tools/hook-scripts/>`_
and put it in

::

    $SVN_DIRECTORY/hooks/

In ``commit-email.pl`` change

::

    my $svnlook = "@SVN_BINDIR@/svnlook";

to

::

    my $svnlook = "/usr/bin/svnlook";

Then create

::

    $SVN_DIRECTORY/hooks/post-commit

which contains (adapt $SVN\_DIRECTORY and email address)

::

    #!/bin/sh

    REPOS="$1"
    REV="$2"

    $SVN_DIRECTORY/hooks/commit-email.pl "$REPOS" "$REV" --from "DIRAC SVN" commit-watchers@example.org

and set

::

    chmod 700 $SVN_DIRECTORY/hooks/post-commit

Retrieved from
"`http://wiki.chem.vu.nl/dirac/index.php/CVS2SVN <http://wiki.chem.vu.nl/dirac/index.php/CVS2SVN>`_"

Views
'''''

-  `Article </dirac/index.php/CVS2SVN>`_
-  `Discussion </dirac/index.php?title=Talk:CVS2SVN&action=edit>`_
-  `View source </dirac/index.php?title=CVS2SVN&action=edit>`_
-  `History </dirac/index.php?title=CVS2SVN&action=history>`_

Personal tools
''''''''''''''

-  `Log in / create
   account </dirac/index.php?title=Special:Userlogin&returnto=CVS2SVN>`_

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

-  `What links here </dirac/index.php/Special:Whatlinkshere/CVS2SVN>`_
-  `Related
   changes </dirac/index.php/Special:Recentchangeslinked/CVS2SVN>`_
-  `Upload file </dirac/index.php/Special:Upload>`_
-  `Special pages </dirac/index.php/Special:Specialpages>`_
-  `Printable version </dirac/index.php?title=CVS2SVN&printable=yes>`_
-  `Permanent link </dirac/index.php?title=CVS2SVN&oldid=7181>`_

`|Powered by MediaWiki| <http://www.mediawiki.org/>`_

-  This page was last modified 07:45, 20 May 2008.
-  This page has been accessed 853 times.
-  `Privacy policy </dirac/index.php/Diracwiki:Privacy_policy>`_
-  `About Diracwiki </dirac/index.php/Diracwiki:About>`_
-  `Disclaimers </dirac/index.php/Diracwiki:General_disclaimer>`_

.. |Powered by
MediaWiki| image:: /dirac/skins/common/images/poweredby_mediawiki_88x31.png
