Manual:WAVE FUNCTION:MVOCAL
===========================

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

Contents
--------

-  `1 \*MVOCAL -- Modified virtual
   orbitals <#.2AMVOCAL_--_Modified_virtual_orbitals>`_

   -  `1.1 .PRINT <#.PRINT>`_
   -  `1.2 .VECMVO <#.VECMVO>`_
   -  `1.3 .INTFLG <#.INTFLG>`_
   -  `1.4 References <#References>`_

\*MVOCAL -- Modified virtual orbitals
=====================================

The virtual canonical Hartree-Fock orbitals are not the optimal choice
for correlated calculations in which the virtual space is truncated; the
lower orbitals tend to be diffuse since they see the *N* + 1 electron
system. One may exploit the freedom of rotation amongst the virtual
orbitals to obtain orbitals better suited for correlated calculations.
The default option in DIRAC is to generate the modified virtual orbitals
by construction of a Fock operator for the system with some electrons
removed, as first proposed by Bauschlicher
\ :sup:``[1] <#_note-bauschlicher:mvo>`_`\ . The orbital string
specified with `.MVOVEC <#.MVOVEC>`_ indicates what occupied orbitals to
remove in the construction of the Hartree-Fock potential.

*Note* that the resulting modified virtual orbitals are no longer
canonical and can *not* be used in MP2 calculations; they can be used in
the `DIRRCI </dirac/index.php/Manual:DIRRCI>`_ and
`RELCCSD </dirac/index.php/Manual:RELCCSD>`_ modules.

.PRINT
------

Print level.

*Default:*

::

    .PRINT
     0

.VECMVO
-------

Read `orbital string </dirac/index.php/Manual:Orbital_strings>`_ of
orbitals to be included in the Hartree-Fock potential when constructing
the cationic Fock operator.

.INTFLG
-------

Specify what two-electron integrals to include during the construction
of the modified virtual orbitals (see
`.INTFLG </dirac/index.php/Manual:HAMILTONIAN#.INTFLG>`_ under
`\*\*HAMILTONIAN </dirac/index.php/Manual:HAMILTONIAN>`_).

*Default:* `.INTFLG </dirac/index.php/Manual:HAMILTONIAN#.INTFLG>`_ from
`\*\*HAMILTONIAN </dirac/index.php/Manual:HAMILTONIAN>`_.

References
----------

#. `↑ <#_ref-bauschlicher:mvo_0>`_ C. W. Bauschlicher, Jr., *The
   construction of modified virtual orbitals (MVO's) which are suited
   for configuration interaction calculations*, J. Chem. Phys. **72**,
   *880* (1980) `electronic
   version <http://dx.doi.org/10.1063/1.439243>`_.

Retrieved from
"`http://wiki.chem.vu.nl/dirac/index.php/Manual:WAVE\_FUNCTION:MVOCAL <http://wiki.chem.vu.nl/dirac/index.php/Manual:WAVE_FUNCTION:MVOCAL>`_"

Views
'''''

-  `Article </dirac/index.php/Manual:WAVE_FUNCTION:MVOCAL>`_
-  `Discussion </dirac/index.php?title=Talk:Manual:WAVE_FUNCTION:MVOCAL&action=edit>`_
-  `View
   source </dirac/index.php?title=Manual:WAVE_FUNCTION:MVOCAL&action=edit>`_
-  `History </dirac/index.php?title=Manual:WAVE_FUNCTION:MVOCAL&action=history>`_

Personal tools
''''''''''''''

-  `Log in / create
   account </dirac/index.php?title=Special:Userlogin&returnto=Manual:WAVE_FUNCTION:MVOCAL>`_

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
   here </dirac/index.php/Special:Whatlinkshere/Manual:WAVE_FUNCTION:MVOCAL>`_
-  `Related
   changes </dirac/index.php/Special:Recentchangeslinked/Manual:WAVE_FUNCTION:MVOCAL>`_
-  `Upload file </dirac/index.php/Special:Upload>`_
-  `Special pages </dirac/index.php/Special:Specialpages>`_
-  `Printable
   version </dirac/index.php?title=Manual:WAVE_FUNCTION:MVOCAL&printable=yes>`_
-  `Permanent
   link </dirac/index.php?title=Manual:WAVE_FUNCTION:MVOCAL&oldid=11421>`_

`|Powered by MediaWiki| <http://www.mediawiki.org/>`_

-  This page was last modified 09:05, 31 July 2008.
-  This page has been accessed 880 times.
-  `Privacy policy </dirac/index.php/Diracwiki:Privacy_policy>`_
-  `About Diracwiki </dirac/index.php/Diracwiki:About>`_
-  `Disclaimers </dirac/index.php/Diracwiki:General_disclaimer>`_

.. |Powered by
MediaWiki| image:: /dirac/skins/common/images/poweredby_mediawiki_88x31.png
