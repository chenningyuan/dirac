Manual:PROPERTIES:LINEAR RESPONSE
=================================

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

Contents
--------

-  `1 General control statements <#General_control_statements>`_

   -  `1.1 .PRINT <#.PRINT>`_

-  `2 Definition of the linear response
   function <#Definition_of_the_linear_response_function>`_

   -  `2.1 .A OPERATOR <#.A_OPERATOR>`_
   -  `2.2 .B OPERATOR <#.B_OPERATOR>`_
   -  `2.3 .OPERATORS <#.OPERATORS>`_
   -  `2.4 .TRIAB <#.TRIAB>`_
   -  `2.5 .B FREQ <#.B_FREQ>`_
   -  `2.6 .IMAGIN <#.IMAGIN>`_
   -  `2.7 .ALLCMB <#.ALLCMB>`_
   -  `2.8 .UNCOUP <#.UNCOUP>`_

-  `3 Control variational parameters <#Control_variational_parameters>`_

   -  `3.1 .OCCUP <#.OCCUP>`_
   -  `3.2 .VIRTUA <#.VIRTUA>`_
   -  `3.3 .SKIPEE <#.SKIPEE>`_
   -  `3.4 .SKIPEP <#.SKIPEP>`_

-  `4 Control reduced equations <#Control_reduced_equations>`_

   -  `4.1 .MAXITR <#.MAXITR>`_
   -  `4.2 .MAXRED <#.MAXRED>`_
   -  `4.3 .THRESH <#.THRESH>`_

-  `5 Control integral contributions <#Control_integral_contributions>`_

   -  `5.1 .INTFLG <#.INTFLG>`_
   -  `5.2 .CNVINT <#.CNVINT>`_
   -  `5.3 .ITRINT <#.ITRINT>`_

-  `6 Control trial vectors <#Control_trial_vectors>`_

   -  `6.1 .REAXVC <#.REAXVC>`_
   -  `6.2 .XLRNRM <#.XLRNRM>`_

-  `7 Advanced/debug flags <#Advanced.2Fdebug_flags>`_

   -  `7.1 .E2CHEK <#.E2CHEK>`_
   -  `7.2 .ONLYSF <#.ONLYSF>`_
   -  `7.3 .ONLYSG <#.ONLYSG>`_
   -  `7.4 .STERNHEIM <#.STERNHEIM>`_
   -  `7.5 .STERNC <#.STERNC>`_
   -  `7.6 .COMPRESSION <#.COMPRESSION>`_
   -  `7.7 .NOPREC <#.NOPREC>`_
   -  `7.8 .RESFAC <#.RESFAC>`_

-  `8 References <#References>`_

Linear Response module written by T. Saue and H. J. Aa. Jensen
\ :sup:``[1] <#_note-saue-jensen:linres>`_`\ .

The general form of a linear response function in the random phase
approximation (RPA) is \ :sup:``[1] <#_note-saue-jensen:linres>`_`\ 
\ :sup:``[2] <#_note-saue:properties>`_`\ 

.. figure:: /dirac/images/math/3/3/a/33ad8ef666cb0f7d49ecd15b16261f82.png
   :align: center
   :alt:  \\langle\\langle A ; B \\rangle\\rangle\_{\\omega} =
   -\\mathbf{E}\_A^{[1]\\dagger} \\bigl( E^{[2]} - \\omega S^{[2]}
   \\bigr)^{-1} \\mathbf{E}\_B^{[1]}

    \\langle\\langle A ; B \\rangle\\rangle\_{\\omega} =
   -\\mathbf{E}\_A^{[1]\\dagger} \\bigl( E^{[2]} - \\omega S^{[2]}
   \\bigr)^{-1} \\mathbf{E}\_B^{[1]}
where |\\mathbf{E}\_A^{[1]}| and |\\mathbf{E}\_B^{[1]}| are property
gradients corresponding to properties A and B, *E*\ :sup:`[2]`\  is the
molecular Hessian, *S*\ :sup:`[2]`\  is the metric, and ω is a
frequency.

The dimension of *E*\ :sup:`[2]`\  does in general not allow an explicit
inversion of the resolvent | \\bigl( E^{[2]} - \\omega S^{[2]}
\\bigr)^{-1} |, so one resorts to iterative techniques.

We write:

.. figure:: /dirac/images/math/e/5/8/e589ae342edd882b5ae789f3aec83ddc.png
   :align: center
   :alt:  \\bigl( E^{[2]} - \\omega S^{[2]} \\bigr)^{-1}
   \\mathbf{E}\_B^{[1]} = \\mathbf{X}

    \\bigl( E^{[2]} - \\omega S^{[2]} \\bigr)^{-1} \\mathbf{E}\_B^{[1]}
   = \\mathbf{X}
By a slight rearrangement we arrive at the linear response equations:

.. figure:: /dirac/images/math/2/1/b/21b88e00d350c78835867c666a7b7015.png
   :align: center
   :alt:  \\bigl( E^{[2]} - \\omega S^{[2]} \\bigr) \\mathbf{X} =
   \\mathbf{E}\_B^{[1]}

    \\bigl( E^{[2]} - \\omega S^{[2]} \\bigr) \\mathbf{X} =
   \\mathbf{E}\_B^{[1]}
The RPA equation is solved by expanding |\\mathbf{X}| in a set of trial
vectors:

.. figure:: /dirac/images/math/0/7/b/07bb04070d9439e9ff637776dcd5b50b.png
   :align: center
   :alt:  \\mathbf{X} = \\sum\_{i=1} \\mathbf{b}\_i a\_i; Y =
   [\\mathbf{b}\_1 \\mathbf{b}\_2 \\cdots \\mathbf{b}\_n]

    \\mathbf{X} = \\sum\_{i=1} \\mathbf{b}\_i a\_i; Y = [\\mathbf{b}\_1
   \\mathbf{b}\_2 \\cdots \\mathbf{b}\_n]
and solving the reduced equation:

.. figure:: /dirac/images/math/e/1/1/e11f9a918dde36f49e11ade054b1fb3e.png
   :align: center
   :alt:  \\bigl( \\tilde{E}^{[2]} - \\omega \\tilde{S}^{[2]} \\bigr)
   \\mathbf{a} = \\tilde{\\mathbf{E}}\_B^{[1]}

    \\bigl( \\tilde{E}^{[2]} - \\omega \\tilde{S}^{[2]} \\bigr)
   \\mathbf{a} = \\tilde{\\mathbf{E}}\_B^{[1]}
in which appear: | \\tilde{E}^{[2]} = Y^\\dagger E^{[2]} Y;
\\tilde{S}^{[2]} = Y^\\dagger S^{[2]} Y; \\tilde{\\mathbf{E}}\_B^{[1]} =
Y^\\dagger \\mathbf{E}\_B^{[1]} |

General control statements
==========================

.PRINT
------

Print level.

*Default:*

::

    .PRINT
     0

Definition of the linear response function
==========================================

.A OPERATOR
-----------

Specification of the A operator (see the `specification of one-electron
operators </dirac/index.php/Manual:One-electron_operators>`_ for
details).

.B OPERATOR
-----------

Specification of the B operator (see the `specification of one-electron
operators </dirac/index.php/Manual:One-electron_operators>`_ for
details).

.OPERATORS
----------

Specification of both the A and B operators (see the `specification of
one-electron
operators </dirac/index.php/Manual:One-electron_operators>`_ for
details).

.TRIAB
------

Enforce triangularity of response function.

.. figure:: /dirac/images/math/2/d/1/2d1ea5d63545e45170ec9f8b4f399a1b.png
   :align: center
   :alt: \\langle\\langle A ; B \\rangle\\rangle\_\\omega \\equiv
   \\langle\\langle B ; A \\rangle\\rangle\_\\omega

   \\langle\\langle A ; B \\rangle\\rangle\_\\omega \\equiv
   \\langle\\langle B ; A \\rangle\\rangle\_\\omega
Only one function is calculated.

*Default:* Deactivated.

.B FREQ
-------

Specify frequencies of operator B.

*Example:* 3 different frequencies.

::

    .B FREQ
     3
     0.001
     0.002
     0.01

*Default:* Static case.

::

    .B FREQ
     1
     0.0

.IMAGIN
-------

Employ imaginary frequencies.

.ALLCMB
-------

Form all possible |\\langle\\langle A ; B \\rangle\\rangle\_\\omega|
even if imaginary.

*Default:* Deactivated.

.UNCOUP
-------

Uncoupled calculation.

Control variational parameters
==============================

.OCCUP
------

For each fermion ircop give an `orbital
string </dirac/index.php/Manual:Orbital_strings>`_ of inactive orbitals
to include in the linear response calculation.

.VIRTUA
-------

For each fermion ircop give an `orbital
string </dirac/index.php/Manual:Orbital_strings>`_ of virtual orbitals
to include in the linear response calculation.

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

Control trial vectors
=====================

.REAXVC
-------

Read solution vectors from file XVCFIL

::

    .REAXVC
    XVCFIL

*Default:* No restart on solution vectors. The file has to have six
characters. Make sure there is no blank character in front of the file
name.

For a restart on solution vectors it is useful to set

::

    .REAXVC
    XVCFIL
    .ITRINT
     0 0

otherwise LS-integrals (and SS-integrals) are switched on later and one
may first iterate away and then back to a possibly converged response
vector.

Often you have a converged SCF wave function along with a response
vector. In this case make sure that

::

    **DIRAC
    #.WAVE FUNCTION

is commented out. Make then also sure that you use the DFCOEF file which
has been obtained in the *same* calculation as the response vector file.
Otherwise you may observe more response solver iterations than
necessary.

.XLRNRM
-------

Normalize trial vectors. Using normalized trial vectors will reduce
efficiency of screening.\ :sup:`[*please
`clarify </dirac/index.php/Clarify>`_*]`\ 

*Default:* Use un-normalized vectors.

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

.STERNHEIM
----------

Set diagonal elements of orbital part of Hessian equal to −
2*m**c*\ :sup:`2`\  for rotations between occupied positive-energy and
virtual negative-energy orbitals.

*Default:* Deactivated.

.STERNC
-------

(Sternheim complement) allows to separate basis set incompleteness from
the replacement of an inner sum over negative-energy orbitals only by
the full sum. In order to benefit from this functionality (only for
specialists !), you should run with print level 2 under properties.

Then you can do a sequence of calculations: 1) .SKIPEP 2) .STERNH 3)
.STERNC The diamagnetic contribution of 1) is the non-relativistic
expectation value, whereas 2) is the Sternheim approximation, that is
replacing orbital energy differences with − 2*m**c*\ :sup:`2`\ . With no
basis set incompleteness the sum of the diamagnetic contribution 2) and
the paramagnetic contribution 3) should equal the diamagnetic
contribution of 1).

*Default:* Deactivated.

.COMPRESSION
------------

Reduce number of orbital variation parameters by checking corresponding
elements of gradient vector against a threshold. This may reduce memory.

*Default:* No compression.

::

    .COMPRESSION
     0.0

.NOPREC
-------

No preconditioning of initial trial vectors.

*Default:* Preconditioning of trial vectors.

.RESFAC
-------

New trial vector will be generated only for variational parameter
classes whose residual has a norm that is larger than a fraction
1/RESFAC of the maximum norm.

*Default:*

::

    .RESFAC
     1000.0

References
==========

#. ↑ \ :sup:``1.0 <#_ref-saue-jensen:linres_0>`_`\ 
   \ :sup:``1.1 <#_ref-saue-jensen:linres_1>`_`\ 
   \ :sup:``1.2 <#_ref-saue-jensen:linres_2>`_`\ 
   \ :sup:``1.3 <#_ref-saue-jensen:linres_3>`_`\  T. Saue and H. J. Aa.
   Jensen, *Linear response at the 4-component relativistic level:
   Application to the frequency-dependent dipole polarizabilities of the
   coinage metal dimers*, J. Chem. Phys. **118**, *522* (2003)
   `electronic version <http://dx.doi.org/10.1063/1.449880>`_.
#. `↑ <#_ref-saue:properties_0>`_ T. Saue, *Post Dirac-Hartree-Fock
   methods - Properties*, in *Relativistic Electronic Structure Theory -
   Part 1: Fundamentals*, edited by P. Schwerdtfeger, Elsevier,
   Amsterdam, 2002.

Retrieved from
"`http://wiki.chem.vu.nl/dirac/index.php/Manual:PROPERTIES:LINEAR\_RESPONSE <http://wiki.chem.vu.nl/dirac/index.php/Manual:PROPERTIES:LINEAR_RESPONSE>`_"

`Category </dirac/index.php/Special:Categories>`_: `Articles needing
clarification </dirac/index.php/Category:Articles_needing_clarification>`_

Views
'''''

-  `Article </dirac/index.php/Manual:PROPERTIES:LINEAR_RESPONSE>`_
-  `Discussion </dirac/index.php/Talk:Manual:PROPERTIES:LINEAR_RESPONSE>`_
-  `View
   source </dirac/index.php?title=Manual:PROPERTIES:LINEAR_RESPONSE&action=edit>`_
-  `History </dirac/index.php?title=Manual:PROPERTIES:LINEAR_RESPONSE&action=history>`_

Personal tools
''''''''''''''

-  `Log in / create
   account </dirac/index.php?title=Special:Userlogin&returnto=Manual:PROPERTIES:LINEAR_RESPONSE>`_

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
   here </dirac/index.php/Special:Whatlinkshere/Manual:PROPERTIES:LINEAR_RESPONSE>`_
-  `Related
   changes </dirac/index.php/Special:Recentchangeslinked/Manual:PROPERTIES:LINEAR_RESPONSE>`_
-  `Upload file </dirac/index.php/Special:Upload>`_
-  `Special pages </dirac/index.php/Special:Specialpages>`_
-  `Printable
   version </dirac/index.php?title=Manual:PROPERTIES:LINEAR_RESPONSE&printable=yes>`_
-  `Permanent
   link </dirac/index.php?title=Manual:PROPERTIES:LINEAR_RESPONSE&oldid=20291>`_

`|Powered by MediaWiki| <http://www.mediawiki.org/>`_

-  This page was last modified 11:40, 1 April 2010.
-  This page has been accessed 1,819 times.
-  `Privacy policy </dirac/index.php/Diracwiki:Privacy_policy>`_
-  `About Diracwiki </dirac/index.php/Diracwiki:About>`_
-  `Disclaimers </dirac/index.php/Diracwiki:General_disclaimer>`_

.. |\\mathbf{E}\_A^{[1]}| image:: /dirac/images/math/6/3/a/63aec7d1787fddef0b15ae2969fc35f5.png
.. |\\mathbf{E}\_B^{[1]}| image:: /dirac/images/math/2/1/8/2182a1c97636f06f80b2642613fe3de6.png
.. | \\bigl( E^{[2]} - \\omega S^{[2]} \\bigr)^{-1}
| image:: /dirac/images/math/6/c/b/6cbb112cbbbb544ae13bedb409423498.png
.. |\\mathbf{X}| image:: /dirac/images/math/5/9/8/598f6444904755dda4a859a1e377468e.png
.. | \\tilde{E}^{[2]} = Y^\\dagger E^{[2]} Y; \\tilde{S}^{[2]} =
Y^\\dagger S^{[2]} Y; \\tilde{\\mathbf{E}}\_B^{[1]} = Y^\\dagger
\\mathbf{E}\_B^{[1]}
| image:: /dirac/images/math/9/b/2/9b2aa6df9f5934136d760e64b1c8fe7d.png
.. |\\langle\\langle A ; B
\\rangle\\rangle\_\\omega| image:: /dirac/images/math/5/d/1/5d1ccb648774ad117a597fb5cb7ad174.png
.. |Powered by
MediaWiki| image:: /dirac/skins/common/images/poweredby_mediawiki_88x31.png
