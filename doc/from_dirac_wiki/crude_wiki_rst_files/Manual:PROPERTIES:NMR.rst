Manual:PROPERTIES:NMR
=====================

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

Contents
--------

-  `1 \*NMR -- Controls evaluation of NMR
   parameters <#.2ANMR_--_Controls_evaluation_of_NMR_parameters>`_

   -  `1.1 .LONDON <#.LONDON>`_
   -  `1.2 .GAUGEO <#.GAUGEO>`_
   -  `1.3 .GO\_ANG <#.GO_ANG>`_
   -  `1.4 .USECM <#.USECM>`_
   -  `1.5 .INTFLG <#.INTFLG>`_
   -  `1.6 .NOTWO <#.NOTWO>`_
   -  `1.7 .NOONEI <#.NOONEI>`_
   -  `1.8 .NOORTH <#.NOORTH>`_
   -  `1.9 .SYMCON <#.SYMCON>`_

\*NMR -- Controls evaluation of NMR parameters
==============================================

This section gives directives for the calculation of NMR parameters.

`|image0| </dirac/index.php/Image:Dirac.png>`_

**Advanced options**

.LONDON
-------

`|image1| </dirac/index.php/Image:Dirac.png>`_

**Only in the development version of DIRAC**

Activate calculations of magnetic properties (NMR shielding constants
and magnetizabilities with London atomic orbitals.

*Default:* Use conventional atomic orbitals.

.GAUGEO
-------

Reads in a user defined gauge origin (in bohr). Note that both the gauge
origin and dipole origin are changed in order to have consistency.

*Default for conventional atomic orbitals:*

::

    .GAUGEO
     0.0 0.0 0.0

.GO\_ANG
--------

Same as `.GAUGEO </dirac/index.php/Manual:PROPERTIES:NMR#.GAUGEO>`_ but
reads coordinates in angstrom.

.USECM
------

Use the center of mass as the gauge origin.

.INTFLG
-------

Specify what two-electron integrals to include in the two-electron
London contributions to the magnetic field property gradient (see
`.INTFLG </dirac/index.php/Manual:HAMILTONIAN#.INTFLG>`_ under
`\*\*HAMILTONIAN </dirac/index.php/Manual:HAMILTONIAN>`_).

*Default:* `.INTFLG </dirac/index.php/Manual:HAMILTONIAN#.INTFLG>`_ from
`\*\*HAMILTONIAN </dirac/index.php/Manual:HAMILTONIAN>`_.

.NOTWO
------

Do not calculate the two-electron London contributions for the magnetic
field property gradient when London atomic orbitals are used.

.NOONEI
-------

Do not calculate the H(0),T(B) reorthonormalization terms for the
magnetic field property gradient when London atomic orbitals are used.

.NOORTH
-------

Do not calculate the T(B),h(mK) reorthonormalization contributions for
the expectation value term when London atomic orbitals are used.

.SYMCON
-------

Employ the symmetric connection for reorthonormalization terms when
using London atomic orbitals.

*Default:* Use the natural connection.

Retrieved from
"`http://wiki.chem.vu.nl/dirac/index.php/Manual:PROPERTIES:NMR <http://wiki.chem.vu.nl/dirac/index.php/Manual:PROPERTIES:NMR>`_"

Views
'''''

-  `Article </dirac/index.php/Manual:PROPERTIES:NMR>`_
-  `Discussion </dirac/index.php?title=Talk:Manual:PROPERTIES:NMR&action=edit>`_
-  `View
   source </dirac/index.php?title=Manual:PROPERTIES:NMR&action=edit>`_
-  `History </dirac/index.php?title=Manual:PROPERTIES:NMR&action=history>`_

Personal tools
''''''''''''''

-  `Log in / create
   account </dirac/index.php?title=Special:Userlogin&returnto=Manual:PROPERTIES:NMR>`_

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
   here </dirac/index.php/Special:Whatlinkshere/Manual:PROPERTIES:NMR>`_
-  `Related
   changes </dirac/index.php/Special:Recentchangeslinked/Manual:PROPERTIES:NMR>`_
-  `Upload file </dirac/index.php/Special:Upload>`_
-  `Special pages </dirac/index.php/Special:Specialpages>`_
-  `Printable
   version </dirac/index.php?title=Manual:PROPERTIES:NMR&printable=yes>`_
-  `Permanent
   link </dirac/index.php?title=Manual:PROPERTIES:NMR&oldid=18321>`_

`|Powered by MediaWiki| <http://www.mediawiki.org/>`_

-  This page was last modified 16:07, 20 September 2009.
-  This page has been accessed 1,371 times.
-  `Privacy policy </dirac/index.php/Diracwiki:Privacy_policy>`_
-  `About Diracwiki </dirac/index.php/Diracwiki:About>`_
-  `Disclaimers </dirac/index.php/Diracwiki:General_disclaimer>`_

.. |image0| image:: /dirac/images/thumb/9/95/Dirac.png/41px-Dirac.png
.. |image1| image:: /dirac/images/thumb/9/95/Dirac.png/41px-Dirac.png
.. |Powered by
MediaWiki| image:: /dirac/skins/common/images/poweredby_mediawiki_88x31.png
