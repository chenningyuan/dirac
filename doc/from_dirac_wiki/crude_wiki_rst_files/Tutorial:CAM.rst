Tutorial:CAM
============

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

In DIRAC, you can calculate with the CAMB3LYP functional using

::

    .DFT
     CAMB3LYP

It is also possible to change the default parameters. The following line
will reproduce the above functional:

::

    .DFT
     CAM p:alpha=0.19 p:beta=0.46 p:mu=0.33 x:slater=1 x:becke=1 c:lyp=0.81 c:vwn5=0.19

but it will give you complete freedom over the parameters. Note that for
this to work correctly you have to apply the following
`patch </dirac/index.php/Known_bugs#User-defined_CAM_functional_line_is_cut_after_character_80>`_.
But as always, please verify your results carefully when changing these
parameters.

Retrieved from
"`http://wiki.chem.vu.nl/dirac/index.php/Tutorial:CAM <http://wiki.chem.vu.nl/dirac/index.php/Tutorial:CAM>`_"

Views
'''''

-  `Article </dirac/index.php/Tutorial:CAM>`_
-  `Discussion </dirac/index.php?title=Talk:Tutorial:CAM&action=edit>`_
-  `View source </dirac/index.php?title=Tutorial:CAM&action=edit>`_
-  `History </dirac/index.php?title=Tutorial:CAM&action=history>`_

Personal tools
''''''''''''''

-  `Log in / create
   account </dirac/index.php?title=Special:Userlogin&returnto=Tutorial:CAM>`_

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
   here </dirac/index.php/Special:Whatlinkshere/Tutorial:CAM>`_
-  `Related
   changes </dirac/index.php/Special:Recentchangeslinked/Tutorial:CAM>`_
-  `Upload file </dirac/index.php/Special:Upload>`_
-  `Special pages </dirac/index.php/Special:Specialpages>`_
-  `Printable
   version </dirac/index.php?title=Tutorial:CAM&printable=yes>`_
-  `Permanent link </dirac/index.php?title=Tutorial:CAM&oldid=26171>`_

`|Powered by MediaWiki| <http://www.mediawiki.org/>`_

-  This page was last modified 21:46, 6 January 2011.
-  This page has been accessed 414 times.
-  `Privacy policy </dirac/index.php/Diracwiki:Privacy_policy>`_
-  `About Diracwiki </dirac/index.php/Diracwiki:About>`_
-  `Disclaimers </dirac/index.php/Diracwiki:General_disclaimer>`_

.. |Powered by
MediaWiki| image:: /dirac/skins/common/images/poweredby_mediawiki_88x31.png
