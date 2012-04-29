Manual:HAMILTONIAN:AMFI
=======================

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

`|image0| </dirac/index.php/Image:Dirac.png>`_

**This page describes unreleased functionality. The keywords may not be
available in your version of DIRAC.**

The AMFI program of B.Schimmelpfennig serves for calculating spin-orbit
contributions to various quasirelativistic two-component Hamiltonians.
Its important part is the independent scalar atomic code - RELSCF (named
as AT34), which is providing atomic orbital coefficients for the AMFI
mean-field summation.

The user may control not only the extent of the output, but also provide
here some parameters for the RELSCF code. Note that closed/open-shell
HF-SCF code gives both the nonrelativistic and the spin-free
second-order Douglass-Krol-Hess SCF energies, which can be compared with
DIRAC results.

RELSCF is using only decontrated basis sets restricted to occupied
shells only; for example, if you set some spdf-basis for a DIRAC
BSS+AMFI run of Magnesium (12 electrons), RELSCF will use, for the sake
of computational effectivity, this basis set with only s- and p-
exponents.

Contents
--------

-  `1 \*AMFI -- AMFI/RELSCF
   directives <#.2AAMFI_--_AMFI.2FRELSCF_directives>`_

   -  `1.1 .PRNT\_A <#.PRNT_A>`_
   -  `1.2 .PRNT\_S <#.PRNT_S>`_
   -  `1.3 .MXITER <#.MXITER>`_
   -  `1.4 .AMFICH <#.AMFICH>`_
   -  `1.5 .NOAMFC <#.NOAMFC>`_
   -  `1.6 .ONLSCF <#.ONLSCF>`_

\*AMFI -- AMFI/RELSCF directives
================================

`|image1| </dirac/index.php/Image:Dirac.png>`_

**Only in the development version of DIRAC**

.PRNT\_A
--------

Print level for AMFI calculations.

*Default:*

::

    .PRNT_A
     0

.PRNT\_S
--------

Print level for the RELSCF part.

*Default:*

::

    .PRNT_S
     0

.MXITER
-------

Maximum number of iterations in RELSCF.

*Default:*

::

    .MXITER
     50

.AMFICH
-------

Set up an artificial charge of the calculated system which has effect
only on the AMFI mean-field summation(s). This does not change the
electronic occupation for the DIRAC SCF procedure.

The reason is that AMFI-attached scalar relativistic Hartree-Fock SCF
module might not be converging for certain high-spin open-shell atoms
(like Pt,Yb), what is anyway expected for a single-determinant code. A
small positive charge (+1, +2) usually helps to get converged molecular
orbitals of given atom with a negligible change in results due to AMFI
mean-field contribution. The AMFI mean-field contributions are large for
core shells, but small for valence orbitals.

In the case of polyatomic systems the overall AMFI-charge is
proportionally distributed over individual atoms and rounded to the
closest integer value.

Example:

::

    .AMFICH
    +1

The user can check the mean-field occupation/charge of each handled atom
in the AMFI/RELSCF output, together with the RELSCF convergence status.

.NOAMFC
-------

Discard selected centers from AMFI calculations. Related the
`.BSS </dirac/index.php/Manual:HAMILTONIAN#.BSS>`_ keyword, whose
numerical value *axyz* must have y=0 (i.e. 'spin-free' from the DIRAC
BSS side, while AMFI provides SO terms of first order in potential).

First line specifies how many atoms are neglected, the second line gives
centers to be omitted. Number for each center of calculated system is
determined in the DIRAC output file.

*Infinite order scalar terms "from the end", and the (one-electron) spin
orbit term (SO1) from AMFI. Centers 1 and 3 are neglected.*

::

    .BSS
    5009
    .NOAMFC
    2
    1  3

.ONLSCF
-------

Stop Dirac calculation after RELSCF program has finished. The
(converged) AO coefficients are written to the file *RELSCF\_COEF* and
can be used to restart from.

This may be useful in combination with the keyword **.AMFICH**, *i.e.*
the user may first run the RELSCF step for an atom with an artificial
charge.

The obtained coefficients can then be used to start the RELSCF step for
the neutral atom for which otherwise one would not obtain convergence.

Input example:

::

    .ONLSCF

Retrieved from
"`http://wiki.chem.vu.nl/dirac/index.php/Manual:HAMILTONIAN:AMFI <http://wiki.chem.vu.nl/dirac/index.php/Manual:HAMILTONIAN:AMFI>`_"

Views
'''''

-  `Article </dirac/index.php/Manual:HAMILTONIAN:AMFI>`_
-  `Discussion </dirac/index.php?title=Talk:Manual:HAMILTONIAN:AMFI&action=edit>`_
-  `View
   source </dirac/index.php?title=Manual:HAMILTONIAN:AMFI&action=edit>`_
-  `History </dirac/index.php?title=Manual:HAMILTONIAN:AMFI&action=history>`_

Personal tools
''''''''''''''

-  `Log in / create
   account </dirac/index.php?title=Special:Userlogin&returnto=Manual:HAMILTONIAN:AMFI>`_

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
   here </dirac/index.php/Special:Whatlinkshere/Manual:HAMILTONIAN:AMFI>`_
-  `Related
   changes </dirac/index.php/Special:Recentchangeslinked/Manual:HAMILTONIAN:AMFI>`_
-  `Upload file </dirac/index.php/Special:Upload>`_
-  `Special pages </dirac/index.php/Special:Specialpages>`_
-  `Printable
   version </dirac/index.php?title=Manual:HAMILTONIAN:AMFI&printable=yes>`_
-  `Permanent
   link </dirac/index.php?title=Manual:HAMILTONIAN:AMFI&oldid=15781>`_

`|Powered by MediaWiki| <http://www.mediawiki.org/>`_

-  This page was last modified 19:48, 10 January 2009.
-  This page has been accessed 1,424 times.
-  `Privacy policy </dirac/index.php/Diracwiki:Privacy_policy>`_
-  `About Diracwiki </dirac/index.php/Diracwiki:About>`_
-  `Disclaimers </dirac/index.php/Diracwiki:General_disclaimer>`_

.. |image0| image:: /dirac/images/thumb/9/95/Dirac.png/41px-Dirac.png
.. |image1| image:: /dirac/images/thumb/9/95/Dirac.png/41px-Dirac.png
.. |Powered by
MediaWiki| image:: /dirac/skins/common/images/poweredby_mediawiki_88x31.png
