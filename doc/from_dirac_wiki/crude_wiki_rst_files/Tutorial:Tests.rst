Tutorial:Tests
==============

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

It is *very important* that you verify that your DIRAC installation
correctly reproduces the reference test set *before* running production
calculations.

To get documentation on the test set, run:

::

    ./runtest --help

As a starter try the quick test that covers most of the standard run
types:

::

    ./runtest --quick

If everything goes well the script will report calculations completed
and provide a summary:

::

     summary:   114 ok

     total time: 00h03m42s

If there is a problem the script will report that one or more tests
failed (wrong output) or crashed (execution error). To learn more go to
the directory **logs** where you will find a text file **\*.txt** with a
short summary of the test that failed. More importantly, there will be a
html-file **\*.htm** to be opened in a browser. Selected parts of your
output are compared to the reference output found in the **result**
subdirectory of your test and a coloring scheme in the html-file will
indicate where things go wrong.

Once you have identified the source of the problem, you may rerun
individual tests with:

::

    ./runtest --tests="test1 test2"

where test1 and test2 should be replaced by the name(s) of the test(s)
that you want to rerun: When everything seems fine try also the complete
test set:

::

    ./runtest --all

Retrieved from
"`http://wiki.chem.vu.nl/dirac/index.php/Tutorial:Tests <http://wiki.chem.vu.nl/dirac/index.php/Tutorial:Tests>`_"

Views
'''''

-  `Article </dirac/index.php/Tutorial:Tests>`_
-  `Discussion </dirac/index.php?title=Talk:Tutorial:Tests&action=edit>`_
-  `View source </dirac/index.php?title=Tutorial:Tests&action=edit>`_
-  `History </dirac/index.php?title=Tutorial:Tests&action=history>`_

Personal tools
''''''''''''''

-  `Log in / create
   account </dirac/index.php?title=Special:Userlogin&returnto=Tutorial:Tests>`_

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
   here </dirac/index.php/Special:Whatlinkshere/Tutorial:Tests>`_
-  `Related
   changes </dirac/index.php/Special:Recentchangeslinked/Tutorial:Tests>`_
-  `Upload file </dirac/index.php/Special:Upload>`_
-  `Special pages </dirac/index.php/Special:Specialpages>`_
-  `Printable
   version </dirac/index.php?title=Tutorial:Tests&printable=yes>`_
-  `Permanent link </dirac/index.php?title=Tutorial:Tests&oldid=31751>`_

`|Powered by MediaWiki| <http://www.mediawiki.org/>`_

-  This page was last modified 19:54, 6 November 2011.
-  This page has been accessed 637 times.
-  `Privacy policy </dirac/index.php/Diracwiki:Privacy_policy>`_
-  `About Diracwiki </dirac/index.php/Diracwiki:About>`_
-  `Disclaimers </dirac/index.php/Diracwiki:General_disclaimer>`_

.. |Powered by
MediaWiki| image:: /dirac/skins/common/images/poweredby_mediawiki_88x31.png
