Manual:PROPERTIES:STEX
======================

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

\*STEX - Core excitation spectra in the Static Exchange approximation
---------------------------------------------------------------------

A well known limitation of the time dependent Hartree-Fock approximation
is that electronic relaxation is not accounted for in the calculation of
excitation energies. This relaxation can be understood as a correlation
effect, since it is due to the simultaneous excitation of an electron
and the relaxation of the remaining *N* − 1 electrons (or rather
orbitals) due to the "hole" left behind by the excited electron. On the
other hand the relaxation is well described already at the Hartree-Fock
level when the excited states are optimized separately, in the so-called
Δ*S**C**F* approach.

The Static Exchange (STEX) approximation is a cheap way to incorporate
relaxation effects into the calculation of core excitation spectra,
where it is particularly important for excitation energies and
transition probabilities. The method works by making a configuration
interaction singles expansion of the excited states with a reference
determinant optimized for core holes in a particular set of orbitals.
The excited states are then not orthogonal to the Hartree-Fock ground
state, and wavefunction overlaps have to be explicitly taken into
account.

A STEX calculation typically proceeds in the following way:

#. The Hartree-Fock groundstate is calculated with Dirac, and the
   wavefunction saved in the ``DFCOEF`` file.
#. The groundstate output is examined and the indices of the core
   orbital of interest are noted.
   If the core orbitals are not well localized to individual atoms it
   may be necessary to add a small
    fractional charge to break a symmetry of the molecule. This change
   can then be removed in the final calculation.
#. An average-of-configuration open shell Hartree-Fock calculation is
   performed to get an approximate wavefunction of the core
    ionized molecule. This wavefunction should then be saved to a file
   ``DFCOEF.ION``.
   In order for the core ionized state to converge it is necessary to
   use overlap selection and possibly orbital freezing.
#. The STEX calculation is run with the specified hole orbitals,
   starting from the previously calculated ``DFCOEF`` and ``DFCOEF.ION``
   files.\ :sup:`[*please provide an
   `example </dirac/index.php/ExampleNeeded>`_*]`\ 

Retrieved from
"`http://wiki.chem.vu.nl/dirac/index.php/Manual:PROPERTIES:STEX <http://wiki.chem.vu.nl/dirac/index.php/Manual:PROPERTIES:STEX>`_"

`Category </dirac/index.php/Special:Categories>`_: `Articles needing
examples </dirac/index.php?title=Category:Articles_needing_examples&action=edit>`_

Views
'''''

-  `Article </dirac/index.php/Manual:PROPERTIES:STEX>`_
-  `Discussion </dirac/index.php?title=Talk:Manual:PROPERTIES:STEX&action=edit>`_
-  `View
   source </dirac/index.php?title=Manual:PROPERTIES:STEX&action=edit>`_
-  `History </dirac/index.php?title=Manual:PROPERTIES:STEX&action=history>`_

Personal tools
''''''''''''''

-  `Log in / create
   account </dirac/index.php?title=Special:Userlogin&returnto=Manual:PROPERTIES:STEX>`_

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
   here </dirac/index.php/Special:Whatlinkshere/Manual:PROPERTIES:STEX>`_
-  `Related
   changes </dirac/index.php/Special:Recentchangeslinked/Manual:PROPERTIES:STEX>`_
-  `Upload file </dirac/index.php/Special:Upload>`_
-  `Special pages </dirac/index.php/Special:Specialpages>`_
-  `Printable
   version </dirac/index.php?title=Manual:PROPERTIES:STEX&printable=yes>`_
-  `Permanent
   link </dirac/index.php?title=Manual:PROPERTIES:STEX&oldid=10161>`_

`|Powered by MediaWiki| <http://www.mediawiki.org/>`_

-  This page was last modified 16:05, 5 June 2008.
-  This page has been accessed 903 times.
-  `Privacy policy </dirac/index.php/Diracwiki:Privacy_policy>`_
-  `About Diracwiki </dirac/index.php/Diracwiki:About>`_
-  `Disclaimers </dirac/index.php/Diracwiki:General_disclaimer>`_

.. |Powered by
MediaWiki| image:: /dirac/skins/common/images/poweredby_mediawiki_88x31.png
