Manual:GOSCIP:Sample inputs
===========================

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

Calculate the spin-orbit components of the 2P ground state of the
fluorine atom using a COSCI wavefunction. This will represent the
"uncorrelated" result since the CI will merely give the energies of the
individual determinants in this simple case. The same calculation can
also be done automatically by specifying
`.RESOLVE </dirac/index.php/Manual:WAVE_FUNCTION:RESOLVE>`_ in the `WAVE
FUNCTION </dirac/index.php/Manual:WAVE_FUNCTION>`_ section of the input.

The full calculation is given as input example 6.cosci.energy in the
test directory:

::

     &GOSCIP NELEC=5, IPRNT=1 &END
     &POPANA THRESH=1.0D-4, DEGEN=1.0D-12, SELPOP=50.0 &END

**NOTE: In the developer's version the first line reads as**

::

     &GOSCIP NELACT=5, IPRNT=1 &END
     &POPANA THRESH=1.0D-4, DEGEN=1.0D-12, SELPOP=50.0 &END

Retrieved from
"`http://wiki.chem.vu.nl/dirac/index.php/Manual:GOSCIP:Sample\_inputs <http://wiki.chem.vu.nl/dirac/index.php/Manual:GOSCIP:Sample_inputs>`_"

Views
'''''

-  `Article </dirac/index.php/Manual:GOSCIP:Sample_inputs>`_
-  `Discussion </dirac/index.php?title=Talk:Manual:GOSCIP:Sample_inputs&action=edit>`_
-  `View
   source </dirac/index.php?title=Manual:GOSCIP:Sample_inputs&action=edit>`_
-  `History </dirac/index.php?title=Manual:GOSCIP:Sample_inputs&action=history>`_

Personal tools
''''''''''''''

-  `Log in / create
   account </dirac/index.php?title=Special:Userlogin&returnto=Manual:GOSCIP:Sample_inputs>`_

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
   here </dirac/index.php/Special:Whatlinkshere/Manual:GOSCIP:Sample_inputs>`_
-  `Related
   changes </dirac/index.php/Special:Recentchangeslinked/Manual:GOSCIP:Sample_inputs>`_
-  `Upload file </dirac/index.php/Special:Upload>`_
-  `Special pages </dirac/index.php/Special:Specialpages>`_
-  `Printable
   version </dirac/index.php?title=Manual:GOSCIP:Sample_inputs&printable=yes>`_
-  `Permanent
   link </dirac/index.php?title=Manual:GOSCIP:Sample_inputs&oldid=19201>`_

`|Powered by MediaWiki| <http://www.mediawiki.org/>`_

-  This page was last modified 13:23, 6 December 2009.
-  This page has been accessed 960 times.
-  `Privacy policy </dirac/index.php/Diracwiki:Privacy_policy>`_
-  `About Diracwiki </dirac/index.php/Diracwiki:About>`_
-  `Disclaimers </dirac/index.php/Diracwiki:General_disclaimer>`_

.. |Powered by
MediaWiki| image:: /dirac/skins/common/images/poweredby_mediawiki_88x31.png
