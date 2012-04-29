Manual:PROPERTIES:EXCITATION ENERGIES
=====================================

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

Calculate excitation energies using time dependent Hartree-Fock or DFT.
The excitation energies are found as the lowest generalized eigenvalues
of the electronic Hessian. DIRAC supports TDDFT kernels from all ground
state functionals included in the code. Currently the iterative
eigenvalue solver may fail to converge more than about twenty roots per
symmetry.

Contents
--------

-  `1 Define excitations and transition
   moments <#Define_excitations_and_transition_moments>`_

   -  `1.1 .EXCITA <#.EXCITA>`_
   -  `1.2 .A OPERATOR <#.A_OPERATOR>`_
   -  `1.3 .ANALYZE <#.ANALYZE>`_

-  `2 Control variational parameters <#Control_variational_parameters>`_

   -  `2.1 .OCCUP <#.OCCUP>`_
   -  `2.2 .VIRTUA <#.VIRTUA>`_
   -  `2.3 .SKIPEE <#.SKIPEE>`_
   -  `2.4 .SKIPEP <#.SKIPEP>`_

-  `3 Control reduced equations <#Control_reduced_equations>`_

   -  `3.1 .MAXITR <#.MAXITR>`_
   -  `3.2 .MAXRED <#.MAXRED>`_
   -  `3.3 .THRESH <#.THRESH>`_

-  `4 Control integral contributions <#Control_integral_contributions>`_

   -  `4.1 .INTFLG <#.INTFLG>`_
   -  `4.2 .CNVINT <#.CNVINT>`_
   -  `4.3 .ITRINT <#.ITRINT>`_

-  `5 Advanced/debug flags <#Advanced.2Fdebug_flags>`_

   -  `5.1 .E2CHEK <#.E2CHEK>`_
   -  `5.2 .ONLYSF <#.ONLYSF>`_
   -  `5.3 .ONLYSG <#.ONLYSG>`_

-  `6 References <#References>`_

Define excitations and transition moments
=========================================

.EXCITA
-------

::

    .EXCITA
    SYM N

Number of excitation energies N calculated in boson symmetry no. SYM.
This keyword can be repeated if you want excitation energies in more
than one boson symmetry.

.A OPERATOR
-----------

Specification of a transition moment operator (see the `specification of
one-electron
operators </dirac/index.php/Manual:One-electron_operators>`_ for
details). This keyword can be given multiple times to add more
operators.

.ANALYZE
--------

Analyze solution vectors and show the most important excitations at the
orbital level.

Control variational parameters
==============================

.OCCUP
------

For each fermion ircop give an `orbital
string </dirac/index.php/Manual:Orbital_strings>`_ of inactive orbitals
from which excitations are allowed. By default excitations from all
occupied orbitals are included in the generalized eigenvalue problem.

Example:

::

    .OCCUP
    1..3
    7,8

This would include excitations from gerade orbitals 1,2,3, and ungerade
orbitals 7 and 8.

.VIRTUA
-------

For each fermion ircop give an `orbital
string </dirac/index.php/Manual:Orbital_strings>`_ of virtual orbitals
to which excitations are allowed. By default excitations to all virtal
orbitals are included in the generalized eigenvalue problem.

.SKIPEE
-------

Exclude all rotations between occupied positive-energy and virtual
positive-energy orbitals.

.SKIPEP
-------

Exclude all rotations between occupied positive-energy and virtual
negative-energy orbitals.

Control reduced equations
=========================

.MAXITR
-------

Maximum number of iterations.

*Default:*

::

    .MAXITR
     30

.MAXRED
-------

Maximum dimension of matrix in reduced system.

*Default:*

::

    .MAXRED
     200

.THRESH
-------

Threshold for convergence of reduced system.

*Default:*

::

    .THRESH
     1.0D-5

Control integral contributions
==============================

The user is encouraged to experiment with these options since they may
have an important effect on run time.

.INTFLG
-------

Specify what two-electron integrals to include (see
`.INTFLG </dirac/index.php/Manual:HAMILTONIAN#.INTFLG>`_ under
`\*\*HAMILTONIAN </dirac/index.php/Manual:HAMILTONIAN>`_).

*Default:* `.INTFLG </dirac/index.php/Manual:HAMILTONIAN#.INTFLG>`_ from
`\*\*HAMILTONIAN </dirac/index.php/Manual:HAMILTONIAN>`_.

.CNVINT
-------

Set threshold for convergence before adding SL and SS integrals to
SCF-iterations.

*2 (real) Arguments:*

::

    .CNVINT
     CNVXQR(1) CNVXQR(2)

*Default:* Very large numbers.

.ITRINT
-------

Set the number of iterations before adding SL and SS integrals to
SCF-iterations.

*Default:*

::

    .ITRINT
     1 1

Advanced/debug flags
====================

.E2CHEK
-------

Generate a complete set of trial vector which implicitly allows the
explicit construction of the electronic Hessian. Only to be used for
small systems !

.ONLYSF
-------

Only call FMOLI in sigmavector routine: only generate one-index
transformed Fock matrix \ :sup:``[1] <#_note-saue-jensen:linres>`_`\ .

.ONLYSG
-------

Only call FMOLI in sigmavector routine: 2-electron Fock matrices using
one-index transformed densities
\ :sup:``[1] <#_note-saue-jensen:linres>`_`\ .

References
==========

#. `↑ <#_ref-saue-jensen:linres_1>`_ **Cite error: Invalid ``<ref>``
   tag; no text was provided for refs named ``saue-jensen:linres``**

Retrieved from
"`http://wiki.chem.vu.nl/dirac/index.php/Manual:PROPERTIES:EXCITATION\_ENERGIES <http://wiki.chem.vu.nl/dirac/index.php/Manual:PROPERTIES:EXCITATION_ENERGIES>`_"

Views
'''''

-  `Article </dirac/index.php/Manual:PROPERTIES:EXCITATION_ENERGIES>`_
-  `Discussion </dirac/index.php?title=Talk:Manual:PROPERTIES:EXCITATION_ENERGIES&action=edit>`_
-  `View
   source </dirac/index.php?title=Manual:PROPERTIES:EXCITATION_ENERGIES&action=edit>`_
-  `History </dirac/index.php?title=Manual:PROPERTIES:EXCITATION_ENERGIES&action=history>`_

Personal tools
''''''''''''''

-  `Log in / create
   account </dirac/index.php?title=Special:Userlogin&returnto=Manual:PROPERTIES:EXCITATION_ENERGIES>`_

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
   here </dirac/index.php/Special:Whatlinkshere/Manual:PROPERTIES:EXCITATION_ENERGIES>`_
-  `Related
   changes </dirac/index.php/Special:Recentchangeslinked/Manual:PROPERTIES:EXCITATION_ENERGIES>`_
-  `Upload file </dirac/index.php/Special:Upload>`_
-  `Special pages </dirac/index.php/Special:Specialpages>`_
-  `Printable
   version </dirac/index.php?title=Manual:PROPERTIES:EXCITATION_ENERGIES&printable=yes>`_
-  `Permanent
   link </dirac/index.php?title=Manual:PROPERTIES:EXCITATION_ENERGIES&oldid=33491>`_

`|Powered by MediaWiki| <http://www.mediawiki.org/>`_

-  This page was last modified 14:24, 27 January 2012.
-  This page has been accessed 806 times.
-  `Privacy policy </dirac/index.php/Diracwiki:Privacy_policy>`_
-  `About Diracwiki </dirac/index.php/Diracwiki:About>`_
-  `Disclaimers </dirac/index.php/Diracwiki:General_disclaimer>`_

.. |Powered by
MediaWiki| image:: /dirac/skins/common/images/poweredby_mediawiki_88x31.png
