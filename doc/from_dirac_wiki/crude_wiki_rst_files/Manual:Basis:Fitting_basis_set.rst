Manual:Basis:Fitting basis set
==============================

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

`|image0| </dirac/index.php/Image:Dirac.png>`_

**This page describes unreleased functionality. The keywords may not be
available in your version of DIRAC.**

Fitting basis set
=================

In the development version of DIRAC it is now possible to make use of
density fitting. In such calculations, it is mandatory to include a
basis set used to expand the density.

This is done adding a line with the FTSET keyword to the MOLECULE file,
after the coordinates for a given atom type

::

    FTSET [options]

Basically all the options available in the definition of large or small
basis sets are applicable to the fit set (i.e. the use of a set from the
library, explicitly typed basis, even/well-tempered etc). Here are some
examples:

Basis sets from the basis set library
-------------------------------------

::

    DIRAC
    Water
    aug-cc-pVDZ basis (note that a tight p has NOT been added)
    C   2    2  X  Y
            8.    1
    O      .0000000000        0.0000000000        -.2249058930
    LARGE BASIS cc-pVDZ
    FTSET BASIS Ahlrichs-Coulomb-Fit

Retrieved from
"`http://wiki.chem.vu.nl/dirac/index.php/Manual:Basis:Fitting\_basis\_set <http://wiki.chem.vu.nl/dirac/index.php/Manual:Basis:Fitting_basis_set>`_"

Views
'''''

-  `Article </dirac/index.php/Manual:Basis:Fitting_basis_set>`_
-  `Discussion </dirac/index.php?title=Talk:Manual:Basis:Fitting_basis_set&action=edit>`_
-  `View
   source </dirac/index.php?title=Manual:Basis:Fitting_basis_set&action=edit>`_
-  `History </dirac/index.php?title=Manual:Basis:Fitting_basis_set&action=history>`_

Personal tools
''''''''''''''

-  `Log in / create
   account </dirac/index.php?title=Special:Userlogin&returnto=Manual:Basis:Fitting_basis_set>`_

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
   here </dirac/index.php/Special:Whatlinkshere/Manual:Basis:Fitting_basis_set>`_
-  `Related
   changes </dirac/index.php/Special:Recentchangeslinked/Manual:Basis:Fitting_basis_set>`_
-  `Upload file </dirac/index.php/Special:Upload>`_
-  `Special pages </dirac/index.php/Special:Specialpages>`_
-  `Printable
   version </dirac/index.php?title=Manual:Basis:Fitting_basis_set&printable=yes>`_
-  `Permanent
   link </dirac/index.php?title=Manual:Basis:Fitting_basis_set&oldid=748>`_

`|Powered by MediaWiki| <http://www.mediawiki.org/>`_

-  This page was last modified 16:07, 12 February 2008.
-  This page has been accessed 1,045 times.
-  `Privacy policy </dirac/index.php/Diracwiki:Privacy_policy>`_
-  `About Diracwiki </dirac/index.php/Diracwiki:About>`_
-  `Disclaimers </dirac/index.php/Diracwiki:General_disclaimer>`_

.. |image0| image:: /dirac/images/thumb/9/95/Dirac.png/41px-Dirac.png
.. |Powered by
MediaWiki| image:: /dirac/skins/common/images/poweredby_mediawiki_88x31.png
