Manual:PROPERTIES:MOLGRD
========================

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

Contents
--------

-  `1 \*MOLGRD -- Controls evaluation of the molecular
   gradient <#.2AMOLGRD_--_Controls_evaluation_of_the_molecular_gradient>`_

   -  `1.1 .TRICK <#.TRICK>`_
   -  `1.2 .INTFLG <#.INTFLG>`_
   -  `1.3 .SCREEN <#.SCREEN>`_
   -  `1.4 .PRINT <#.PRINT>`_
   -  `1.5 .NUMGRA <#.NUMGRA>`_

\*MOLGRD -- Controls evaluation of the molecular gradient
=========================================================

The section gives directives for control of a single point evaluation of
the molecular gradient activated with
`.MOLGRD </dirac/index.php/Manual:PROPERTIES#.MOLGRD>`_ under
`\*\*PROPERTIES </dirac/index.php/Manual:PROPERTIES>`_.

The molecular gradient module is usually called from within the geometry
optimization module.

`|image0| </dirac/index.php/Image:Dirac.png>`_

**Advanced options**

.TRICK
------

If activated the calculation of the SL and SS two-electron integral
contributions to the gradient is skipped if their contribution is
estimated to be small. An "empirical" estimate of the norms of the LS
and SS two-electron gradients based on the norm of the LL two-electron
gradient. If deactivated all contributions to gradient are calculated
and if activated then only "necessary" contributions are calculated.

*Default:* Deactivated.

.INTFLG
-------

Specify what two-electron integrals to include (see
`.INTFLG </dirac/index.php/Manual:HAMILTONIAN#.INTFLG>`_ under
`\*\*HAMILTONIAN </dirac/index.php/Manual:HAMILTONIAN>`_).

*Default:* `.INTFLG </dirac/index.php/Manual:HAMILTONIAN#.INTFLG>`_ from
`\*\*HAMILTONIAN </dirac/index.php/Manual:HAMILTONIAN>`_.

.SCREEN
-------

Screening threshold.

*Default:*

::

    .SCREEN
     1.0D-17

Which means that screening is turned off by default. This is due to an
enormous demand of memory.

`|image1| </dirac/index.php/Image:HJAaJ.png>`_

**Programmers options**

.PRINT
------

Print level.

A print level of less than 3 gives only the total gradient. Print levels
from 3+ gives individual contributions to the gradient (kinetic energy
gradient, nuclear attraction gradient etc.). Print levels of 5+ print
some matrices, and 10+ give massive output!

*Default:* Print level from
`\*\*GENERAL </dirac/index.php/Manual:GENERAL>`_.

.NUMGRA
-------

Calculate the molecular gradient using numerical differentiation.

*Default:* Analytical evaluation if implemented.

Retrieved from
"`http://wiki.chem.vu.nl/dirac/index.php/Manual:PROPERTIES:MOLGRD <http://wiki.chem.vu.nl/dirac/index.php/Manual:PROPERTIES:MOLGRD>`_"

Views
'''''

-  `Article </dirac/index.php/Manual:PROPERTIES:MOLGRD>`_
-  `Discussion </dirac/index.php?title=Talk:Manual:PROPERTIES:MOLGRD&action=edit>`_
-  `View
   source </dirac/index.php?title=Manual:PROPERTIES:MOLGRD&action=edit>`_
-  `History </dirac/index.php?title=Manual:PROPERTIES:MOLGRD&action=history>`_

Personal tools
''''''''''''''

-  `Log in / create
   account </dirac/index.php?title=Special:Userlogin&returnto=Manual:PROPERTIES:MOLGRD>`_

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
   here </dirac/index.php/Special:Whatlinkshere/Manual:PROPERTIES:MOLGRD>`_
-  `Related
   changes </dirac/index.php/Special:Recentchangeslinked/Manual:PROPERTIES:MOLGRD>`_
-  `Upload file </dirac/index.php/Special:Upload>`_
-  `Special pages </dirac/index.php/Special:Specialpages>`_
-  `Printable
   version </dirac/index.php?title=Manual:PROPERTIES:MOLGRD&printable=yes>`_
-  `Permanent
   link </dirac/index.php?title=Manual:PROPERTIES:MOLGRD&oldid=3841>`_

`|Powered by MediaWiki| <http://www.mediawiki.org/>`_

-  This page was last modified 12:12, 10 March 2008.
-  This page has been accessed 821 times.
-  `Privacy policy </dirac/index.php/Diracwiki:Privacy_policy>`_
-  `About Diracwiki </dirac/index.php/Diracwiki:About>`_
-  `Disclaimers </dirac/index.php/Diracwiki:General_disclaimer>`_

.. |image0| image:: /dirac/images/thumb/9/95/Dirac.png/41px-Dirac.png
.. |image1| image:: /dirac/images/thumb/6/65/HJAaJ.png/41px-HJAaJ.png
.. |Powered by
MediaWiki| image:: /dirac/skins/common/images/poweredby_mediawiki_88x31.png
