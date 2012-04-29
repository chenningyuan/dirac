Manual:KR-CI properties
=======================

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

*2- and 4-component relativistic KR-CI property module* written by
Stefan Knecht and Hans Joergen Aa. Jensen, parallelization by Stefan
Knecht

The KR-CI property module \ :sup:``[1] <#_note-knecht:thesis>`_`\  takes
advantage of the 2- and 4-component KR-CI module. It can be used for the
computation of:

::

    - permanent dipole moments in ground and excited states, 

::

    - transition dipole moments 

::

    - computation of ,  and  (Ω values) expectation values (only for linear molecules).

Upon request the analysis of other one-electron properties may also be
implemented. The module could in principle be used for each one-electron
operator that is specified in
`[1] <http://wiki.chem.vu.nl/dirac/index.php/Manual:One-electron_operators>`_.

**NOTE: the following keywords must be placed under the input deck**

::

    *KRCICALC

together with the appropriate CI calculation input (GASSH, GASSPC, etc.)

Contents
--------

-  `1 implemented keywords <#implemented_keywords>`_

   -  `1.1 .DIPMOM <#.DIPMOM>`_
   -  `1.2 .TRDM <#.TRDM>`_
   -  `1.3 .OMEGAQ <#.OMEGAQ>`_

-  `2 References <#References>`_

implemented keywords
====================

.DIPMOM
-------

compute the permanent dipole moments in electronic ground and excited
states.

.TRDM
-----

compute the transition dipole moments between electronic states.

.OMEGAQ
-------

compute the expectation values of the spin- and angular momentum
operator in z-direction and print the total |<\\hat{j}\_z> =
<\\hat{s}\_z> + <\\hat{l}\_z>| expectation value (**Ω** value) for each
electronic state.

::

    Note: this applies only to linear molecules.

References
==========

#. `↑ <#_ref-knecht:thesis_0>`_ S. R. Knecht, *PhD thesis*, University
   of Duesseldorf (Germany) '***, '** (2009)* .

Retrieved from
"`http://wiki.chem.vu.nl/dirac/index.php/Manual:KR-CI\_properties <http://wiki.chem.vu.nl/dirac/index.php/Manual:KR-CI_properties>`_"

Views
'''''

-  `Article </dirac/index.php/Manual:KR-CI_properties>`_
-  `Discussion </dirac/index.php?title=Talk:Manual:KR-CI_properties&action=edit>`_
-  `View
   source </dirac/index.php?title=Manual:KR-CI_properties&action=edit>`_
-  `History </dirac/index.php?title=Manual:KR-CI_properties&action=history>`_

Personal tools
''''''''''''''

-  `Log in / create
   account </dirac/index.php?title=Special:Userlogin&returnto=Manual:KR-CI_properties>`_

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
   here </dirac/index.php/Special:Whatlinkshere/Manual:KR-CI_properties>`_
-  `Related
   changes </dirac/index.php/Special:Recentchangeslinked/Manual:KR-CI_properties>`_
-  `Upload file </dirac/index.php/Special:Upload>`_
-  `Special pages </dirac/index.php/Special:Specialpages>`_
-  `Printable
   version </dirac/index.php?title=Manual:KR-CI_properties&printable=yes>`_
-  `Permanent
   link </dirac/index.php?title=Manual:KR-CI_properties&oldid=24061>`_

`|Powered by MediaWiki| <http://www.mediawiki.org/>`_

-  This page was last modified 17:36, 11 October 2010.
-  This page has been accessed 425 times.
-  `Privacy policy </dirac/index.php/Diracwiki:Privacy_policy>`_
-  `About Diracwiki </dirac/index.php/Diracwiki:About>`_
-  `Disclaimers </dirac/index.php/Diracwiki:General_disclaimer>`_

.. |<\\hat{j}\_z> = <\\hat{s}\_z> +
<\\hat{l}\_z>| image:: /dirac/images/math/6/6/d/66d28911e2198a24014da6b252f5fc93.png
.. |Powered by
MediaWiki| image:: /dirac/skins/common/images/poweredby_mediawiki_88x31.png
