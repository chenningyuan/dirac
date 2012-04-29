Manual:INTEGRALS:READIN
=======================

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

Contents
--------

-  `1 \*READIN -- Reading basis set
   file <#.2AREADIN_--_Reading_basis_set_file>`_

   -  `1.1 .UNCONTRACT <#.UNCONTRACT>`_
   -  `1.2 .PRINT <#.PRINT>`_
   -  `1.3 .MAXPRI <#.MAXPRI>`_

\*READIN -- Reading basis set file
==================================

This subsection allows changes of defaults in the reading of the `basis
file </dirac/index.php/Manual:Basis:Title_and_symmetry>`_.

`|image0| </dirac/index.php/Image:Dirac.png>`_

**Advanced options**

.UNCONTRACT
-----------

Decontract basis sets specified as contracted (when being read from the
library, for example). In the case of using non-relativistically
contracted basis set the decontraction is necessary for heavy elements.
Two-component quasirelativistic Hamiltonians (like IOTC) work only with
decontracted basis set.

*Default:* Decontract only the small component.

`|image1| </dirac/index.php/Image:HJAaJ.png>`_

**Programmers options**

.PRINT
------

Print level in the reading of the basis file.

*Default:* Print level from
`\*\*INTEGRALS </dirac/index.php/Manual:INTEGRALS>`_.

.MAXPRI
-------

Maximum number of primitive functions in a given block in basis file.

*Default:*

::

    .MAXPRI
     22

Retrieved from
"`http://wiki.chem.vu.nl/dirac/index.php/Manual:INTEGRALS:READIN <http://wiki.chem.vu.nl/dirac/index.php/Manual:INTEGRALS:READIN>`_"

Views
'''''

-  `Article </dirac/index.php/Manual:INTEGRALS:READIN>`_
-  `Discussion </dirac/index.php?title=Talk:Manual:INTEGRALS:READIN&action=edit>`_
-  `View
   source </dirac/index.php?title=Manual:INTEGRALS:READIN&action=edit>`_
-  `History </dirac/index.php?title=Manual:INTEGRALS:READIN&action=history>`_

Personal tools
''''''''''''''

-  `Log in / create
   account </dirac/index.php?title=Special:Userlogin&returnto=Manual:INTEGRALS:READIN>`_

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
   here </dirac/index.php/Special:Whatlinkshere/Manual:INTEGRALS:READIN>`_
-  `Related
   changes </dirac/index.php/Special:Recentchangeslinked/Manual:INTEGRALS:READIN>`_
-  `Upload file </dirac/index.php/Special:Upload>`_
-  `Special pages </dirac/index.php/Special:Specialpages>`_
-  `Printable
   version </dirac/index.php?title=Manual:INTEGRALS:READIN&printable=yes>`_
-  `Permanent
   link </dirac/index.php?title=Manual:INTEGRALS:READIN&oldid=11011>`_

`|Powered by MediaWiki| <http://www.mediawiki.org/>`_

-  This page was last modified 15:50, 12 June 2008.
-  This page has been accessed 1,233 times.
-  `Privacy policy </dirac/index.php/Diracwiki:Privacy_policy>`_
-  `About Diracwiki </dirac/index.php/Diracwiki:About>`_
-  `Disclaimers </dirac/index.php/Diracwiki:General_disclaimer>`_

.. |image0| image:: /dirac/images/thumb/9/95/Dirac.png/41px-Dirac.png
.. |image1| image:: /dirac/images/thumb/6/65/HJAaJ.png/41px-HJAaJ.png
.. |Powered by
MediaWiki| image:: /dirac/skins/common/images/poweredby_mediawiki_88x31.png
