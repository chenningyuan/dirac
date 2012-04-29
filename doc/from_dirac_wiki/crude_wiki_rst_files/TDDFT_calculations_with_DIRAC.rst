TDDFT calculations with DIRAC
=============================

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

Here I would like to give a typical input for a TDDFT calculation of
excitation energies with DIRAC and discuss the important keywords
(colored).

This is a TDDFT calculation of the Zn atom with the PBE functional.

To obtain good TDDFT excitation energies it is *essential* to include
the spin density contribution (spin magnetization contribution). This is
done with the .SDFT keyword. Without spin density, triplet excitation
energies correspond to KS orbital energy differences.

We suppress orbital rotations between positive and negative energy
orbitals with .SKIPEP.

Also include the .TRIPLET keyword. This keyword is without effect (but
also without harm) when spin-orbit coupling is turned on. However,
without spin-orbit coupling (nonrelativistic or spin-free calculations)
this keyword allows the spin density contribution in spin-forbidden
transitions. Without spin-orbit coupling and without .TRIPLET, triplet
energies again correspond to KS orbital energy differences. This keyword
does *not* mean that only triplet states are solved for.

::

    **DIRAC
    .WAVE FUNCTION
    .PROPERTIES
    .ANALYZE
    **HAMILTONIAN
    .LVCORR
    .DFT
     PBE
    *DFT
    .SDFT
    **INTEGRALS
    *READIN
    .UNCONTRACT
    **WAVE FUNCTION
    .DHF
    *DHFCAL
    .EVCCNV
     1.0D-9 1.0D-6
    .CLOSED
     18 12      !zinc atom
    **PROPERTIES
    *EXCITATION ENERGIES
    .TRIPLET
    .SKIPEP
    .ANALYZE
    .EXCITATION ENERGIES
     3 4
    .EXCITATION ENERGIES
     5 4
    .EXCITATION ENERGIES
     8 4
    .OPERATOR
     XDIPLEN
    .OPERATOR
     YDIPLEN
    .OPERATOR
     ZDIPLEN
    **ANALYZE
    .MULPOP
    *MULPOP
    .VECPOP
     1..oo
     1..oo      !remove this line if there is only one Fermion irreducible corepresentation
    *END OF

With

::

    .EXCITATION ENERGIES
     3 4
    .EXCITATION ENERGIES
     5 4
    .EXCITATION ENERGIES
     8 4

we ask for the four lowest excitation energies in the irreps 3, 5, and
8.

In this example (Zn, *D*\ :sub:`2*h*`\ ) these are the ungerade irreps
B2u, B1u, and Au.

Finally, it is a good idea to always include

::

    .OPERATOR
     XDIPLEN
    .OPERATOR
     YDIPLEN
    .OPERATOR
     ZDIPLEN

to get transition dipole moments and oscillator strengths (if symmetry
allowed).

Retrieved from
"`http://wiki.chem.vu.nl/dirac/index.php/TDDFT\_calculations\_with\_DIRAC <http://wiki.chem.vu.nl/dirac/index.php/TDDFT_calculations_with_DIRAC>`_"

Views
'''''

-  `Article </dirac/index.php/TDDFT_calculations_with_DIRAC>`_
-  `Discussion </dirac/index.php?title=Talk:TDDFT_calculations_with_DIRAC&action=edit>`_
-  `View
   source </dirac/index.php?title=TDDFT_calculations_with_DIRAC&action=edit>`_
-  `History </dirac/index.php?title=TDDFT_calculations_with_DIRAC&action=history>`_

Personal tools
''''''''''''''

-  `Log in / create
   account </dirac/index.php?title=Special:Userlogin&returnto=TDDFT_calculations_with_DIRAC>`_

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
   here </dirac/index.php/Special:Whatlinkshere/TDDFT_calculations_with_DIRAC>`_
-  `Related
   changes </dirac/index.php/Special:Recentchangeslinked/TDDFT_calculations_with_DIRAC>`_
-  `Upload file </dirac/index.php/Special:Upload>`_
-  `Special pages </dirac/index.php/Special:Specialpages>`_
-  `Printable
   version </dirac/index.php?title=TDDFT_calculations_with_DIRAC&printable=yes>`_
-  `Permanent
   link </dirac/index.php?title=TDDFT_calculations_with_DIRAC&oldid=15381>`_

`|Powered by MediaWiki| <http://www.mediawiki.org/>`_

-  This page was last modified 11:35, 6 December 2008.
-  This page has been accessed 847 times.
-  `Privacy policy </dirac/index.php/Diracwiki:Privacy_policy>`_
-  `About Diracwiki </dirac/index.php/Diracwiki:About>`_
-  `Disclaimers </dirac/index.php/Diracwiki:General_disclaimer>`_

.. |Powered by
MediaWiki| image:: /dirac/skins/common/images/poweredby_mediawiki_88x31.png
