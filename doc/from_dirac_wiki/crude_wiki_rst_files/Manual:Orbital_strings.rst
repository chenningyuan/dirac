Manual:Orbital strings
======================

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

Specification of orbital strings
================================

A subset of orbitals can be specified in two ways. Usually it suffices
to specify the keyword energy, followed by three real numbers to specify
a lower limit, upper limit and quasi-degeneracy threshold. All orbitals
that have an energy between the two thresholds are then selected. The
quasi-degeneracy threshold is used to extend the range if the lower or
upper threshold cuts through a set of closely energy-spaced orbitals. It
will move the upper or lower limit until it encounters a large enough
gap between orbital energies. The default value used in the AO-MO
integral transformation program to specify the active space for
correlated calculations is e.g. given by the string

::

    energy -10.0 20.0 1.0

which picks out all orbitals that have an energy between -10.0 and +20.0
hartree (a.u.), with a minimum gap of 1.0 hartree.

More detailed control is possible by using an orbital string. Here a
character string is given where individual orbitals are listed separated
by a comma or as a range of orbitals. The range limit is then separated
by "..". Positive energy and negative energy orbitals are indicated
using positive and negative numbers, respectively. For instance the
orbital string

::

    -5..5,7

picks out the five upper negative energy ("positronic") and five lower
positive energy ("electronic") orbitals plus positive energy orbital
number seven.

If all orbitals (both positive and negative energy) are wanted, then one
may write

::

    all

The concept of orbitals strings gives great flexibility in the
specification of orbitals and serves as an alternative to editing of
vector files. A drawback is, however, that one needs to count the number
of orbitals manually, which may become cumbersome for cases with large
ranges of virtual orbitals.

Retrieved from
"`http://wiki.chem.vu.nl/dirac/index.php/Manual:Orbital\_strings <http://wiki.chem.vu.nl/dirac/index.php/Manual:Orbital_strings>`_"

Views
'''''

-  `Article </dirac/index.php/Manual:Orbital_strings>`_
-  `Discussion </dirac/index.php?title=Talk:Manual:Orbital_strings&action=edit>`_
-  `View
   source </dirac/index.php?title=Manual:Orbital_strings&action=edit>`_
-  `History </dirac/index.php?title=Manual:Orbital_strings&action=history>`_

Personal tools
''''''''''''''

-  `Log in / create
   account </dirac/index.php?title=Special:Userlogin&returnto=Manual:Orbital_strings>`_

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
   here </dirac/index.php/Special:Whatlinkshere/Manual:Orbital_strings>`_
-  `Related
   changes </dirac/index.php/Special:Recentchangeslinked/Manual:Orbital_strings>`_
-  `Upload file </dirac/index.php/Special:Upload>`_
-  `Special pages </dirac/index.php/Special:Specialpages>`_
-  `Printable
   version </dirac/index.php?title=Manual:Orbital_strings&printable=yes>`_
-  `Permanent
   link </dirac/index.php?title=Manual:Orbital_strings&oldid=22501>`_

`|Powered by MediaWiki| <http://www.mediawiki.org/>`_

-  This page was last modified 10:33, 6 September 2010.
-  This page has been accessed 1,456 times.
-  `Privacy policy </dirac/index.php/Diracwiki:Privacy_policy>`_
-  `About Diracwiki </dirac/index.php/Diracwiki:About>`_
-  `Disclaimers </dirac/index.php/Diracwiki:General_disclaimer>`_

.. |Powered by
MediaWiki| image:: /dirac/skins/common/images/poweredby_mediawiki_88x31.png
