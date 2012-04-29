SVN and Git
===========

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

Motivation
----------

We use the version control system
`Subversion <http://subversion.tigris.org/>`_ (SVN) for DIRAC and there
are good reasons for it.

The idea of this how-to is to show that it is useful and easy to use SVN
in combination with `Git <http://git.or.cz/>`_ for the programming and
testing of DIRAC.

It can be useful if you want to:

-  Have full advantage of a version control system: commit often,
   experiment, revert, ..., and you don't want to worry about breaking
   the main trunk and you don't want to run the test suite before each
   commit.
-  Track updates on the main SVN server but commit your changes to the
   main trunk only from time to time when you feel like it (when you
   have ran the test suite).
-  Test your changes for different builds, on different machines
   **before** committing them to the main trunk.

Why Git? Because Git it makes it very easy to replicate repositories,
create experimental branches, push/pull from/to local or remote
repositories.

Installation
------------

On debian based systems:

::

    sudo apt-get install git-core git-svn

Git repository for DIRAC in two steps
-------------------------------------

git-svn is a very nice tool (see `An introduction to git-svn for
Subversion/SVK users and
deserters <http://utsl.gen.nz/talks/git-svn/intro.html#local-commit>`_)
to interface SVN repositories with Git. It will allow you to track the
main SVN repository and allow to push a series of commits to SVN.

For my development I decided to use a more primitive approach to set up
a Git repository for DIRAC which is described below.

I prefer to update from the main SVN repository "by hand". And I will
not be able to push a series of Git commits to SVN automatically but
again prefer to commit changes to SVN in larger chunks and select them
"by hand".

First set the necessary variables and check out a fresh working copy
from the $SVN\_REPOSITORY to the directory $GIT\_REPOSITORY

::

    svn co $SVN_REPOSITORY $GIT_REPOSITORY

Then

::

    cd $GIT_REPOSITORY

    git init
    git add .
    git commit -m 'initial repository layout'

Done.

Now you can program, have fun, branch, clone, ... use all the nice
tricks Git offers (see the `Git crash course for SVN
users <http://git.or.cz/course/svn.html>`_).

radovan: under development ...

Retrieved from
"`http://wiki.chem.vu.nl/dirac/index.php/SVN\_and\_Git <http://wiki.chem.vu.nl/dirac/index.php/SVN_and_Git>`_"

Views
'''''

-  `Article </dirac/index.php/SVN_and_Git>`_
-  `Discussion </dirac/index.php?title=Talk:SVN_and_Git&action=edit>`_
-  `View source </dirac/index.php?title=SVN_and_Git&action=edit>`_
-  `History </dirac/index.php?title=SVN_and_Git&action=history>`_

Personal tools
''''''''''''''

-  `Log in / create
   account </dirac/index.php?title=Special:Userlogin&returnto=SVN_and_Git>`_

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
   here </dirac/index.php/Special:Whatlinkshere/SVN_and_Git>`_
-  `Related
   changes </dirac/index.php/Special:Recentchangeslinked/SVN_and_Git>`_
-  `Upload file </dirac/index.php/Special:Upload>`_
-  `Special pages </dirac/index.php/Special:Specialpages>`_
-  `Printable
   version </dirac/index.php?title=SVN_and_Git&printable=yes>`_
-  `Permanent link </dirac/index.php?title=SVN_and_Git&oldid=8821>`_

`|Powered by MediaWiki| <http://www.mediawiki.org/>`_

-  This page was last modified 12:17, 29 May 2008.
-  This page has been accessed 160 times.
-  `Privacy policy </dirac/index.php/Diracwiki:Privacy_policy>`_
-  `About Diracwiki </dirac/index.php/Diracwiki:About>`_
-  `Disclaimers </dirac/index.php/Diracwiki:General_disclaimer>`_

.. |Powered by
MediaWiki| image:: /dirac/skins/common/images/poweredby_mediawiki_88x31.png
