Manual:PROPERTIES:QUADRATIC RESPONSE
====================================

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

Contents
--------

-  `1 General control statements <#General_control_statements>`_

   -  `1.1 .PRINT <#.PRINT>`_

-  `2 Definition of the quadratic response
   function <#Definition_of_the_quadratic_response_function>`_

   -  `2.1 .DIPLEN <#.DIPLEN>`_
   -  `2.2 .A OPERATOR <#.A_OPERATOR>`_
   -  `2.3 .B OPERATOR <#.B_OPERATOR>`_
   -  `2.4 .C OPERATOR <#.C_OPERATOR>`_
   -  `2.5 .B FREQ <#.B_FREQ>`_
   -  `2.6 .C FREQ <#.C_FREQ>`_
   -  `2.7 .ALLCMB <#.ALLCMB>`_

-  `3 Excited state properties <#Excited_state_properties>`_

   -  `3.1 .EXCPRP <#.EXCPRP>`_

-  `4 Control variational parameters <#Control_variational_parameters>`_

   -  `4.1 .SKIPEE <#.SKIPEE>`_
   -  `4.2 .SKIPEP <#.SKIPEP>`_

-  `5 Control reduced equations <#Control_reduced_equations>`_

   -  `5.1 .MAXITR <#.MAXITR>`_
   -  `5.2 .MAXRED <#.MAXRED>`_
   -  `5.3 .THRESH <#.THRESH>`_

-  `6 Control integral contributions <#Control_integral_contributions>`_

   -  `6.1 .INTFLG <#.INTFLG>`_
   -  `6.2 .CNVINT <#.CNVINT>`_
   -  `6.3 .ITRINT <#.ITRINT>`_

-  `7 Control trial vectors <#Control_trial_vectors>`_

   -  `7.1 .XQRNRM <#.XQRNRM>`_

-  `8 Advanced/debug flags <#Advanced.2Fdebug_flags>`_

   -  `8.1 .NOPREC <#.NOPREC>`_
   -  `8.2 .RESFAC <#.RESFAC>`_

-  `9 References <#References>`_

This section gives directives for the calculation of quadratic response
functions \ :sup:``[1] <#_note-saue:properties>`_`\  |\\langle\\langle
A ; B , C \\rangle\\rangle\_{\\omega\_1,\\omega\_2}|

General control statements
==========================

.PRINT
------

Print level.

*Default:*

::

    .PRINT
     0

Definition of the quadratic response function
=============================================

.DIPLEN
-------

Specification of dipole operators for A, B, and C (see the
`specification of one-electron
operators </dirac/index.php/Manual:One-electron_operators>`_ for
details).

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

.C OPERATOR
-----------

Specification of the C operator (see the `specification of one-electron
operators </dirac/index.php/Manual:One-electron_operators>`_ for
details).

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

.C FREQ
-------

Specify frequencies of operator C (see `.B FREQ <#.B_FREQ>`_).

.ALLCMB
-------

Evaluate all nonzero quadratic response functions and thereby
disregarding analysis of overall permutational symmetry.

*Default:* Evaluate only unique, nonzero, response functions.

Excited state properties
========================

`|image1| </dirac/index.php/Image:Dirac.png>`_

**This page describes unreleased functionality. The keywords may not be
available in your version of DIRAC.**

First order properties of excited states can be computed from the
quadratic response function.

.EXCPRP
-------

Give the number of "left" and "right" states in each boson symmetry.

*Example*:

::

    .EXCPRP
    5 5 5 5
    0 0 0 0

Compute the transition moments |\\langle i\|\\hat{A}\|i\\rangle|, where
i goes from 1 to 5 in each symmetry (four symmetries in this case). The
zeros can be substituted for positive integers to generate elements
|\\langle i\|\\hat{A}\|j\\rangle|.

Control variational parameters
==============================

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
     100

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

.XQRNRM
-------

Normalize trial vectors. Using normalized trial vectors will reduce
efficiency of screening.

*Default:* Use un-normalized vectors.

Advanced/debug flags
====================

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

#. `↑ <#_ref-saue:properties_0>`_ T. Saue, *Post Dirac-Hartree-Fock
   methods - Properties*, in *Relativistic Electronic Structure Theory -
   Part 1: Fundamentals*, edited by P. Schwerdtfeger, Elsevier,
   Amsterdam, 2002.

Retrieved from
"`http://wiki.chem.vu.nl/dirac/index.php/Manual:PROPERTIES:QUADRATIC\_RESPONSE <http://wiki.chem.vu.nl/dirac/index.php/Manual:PROPERTIES:QUADRATIC_RESPONSE>`_"

Views
'''''

-  `Article </dirac/index.php/Manual:PROPERTIES:QUADRATIC_RESPONSE>`_
-  `Discussion </dirac/index.php?title=Talk:Manual:PROPERTIES:QUADRATIC_RESPONSE&action=edit>`_
-  `View
   source </dirac/index.php?title=Manual:PROPERTIES:QUADRATIC_RESPONSE&action=edit>`_
-  `History </dirac/index.php?title=Manual:PROPERTIES:QUADRATIC_RESPONSE&action=history>`_

Personal tools
''''''''''''''

-  `Log in / create
   account </dirac/index.php?title=Special:Userlogin&returnto=Manual:PROPERTIES:QUADRATIC_RESPONSE>`_

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
   here </dirac/index.php/Special:Whatlinkshere/Manual:PROPERTIES:QUADRATIC_RESPONSE>`_
-  `Related
   changes </dirac/index.php/Special:Recentchangeslinked/Manual:PROPERTIES:QUADRATIC_RESPONSE>`_
-  `Upload file </dirac/index.php/Special:Upload>`_
-  `Special pages </dirac/index.php/Special:Specialpages>`_
-  `Printable
   version </dirac/index.php?title=Manual:PROPERTIES:QUADRATIC_RESPONSE&printable=yes>`_
-  `Permanent
   link </dirac/index.php?title=Manual:PROPERTIES:QUADRATIC_RESPONSE&oldid=19781>`_

`|Powered by MediaWiki| <http://www.mediawiki.org/>`_

-  This page was last modified 15:07, 2 February 2010.
-  This page has been accessed 1,029 times.
-  `Privacy policy </dirac/index.php/Diracwiki:Privacy_policy>`_
-  `About Diracwiki </dirac/index.php/Diracwiki:About>`_
-  `Disclaimers </dirac/index.php/Diracwiki:General_disclaimer>`_

.. |\\langle\\langle A ; B , C
\\rangle\\rangle\_{\\omega\_1,\\omega\_2}| image:: /dirac/images/math/e/b/0/eb02ef7a48d65d5a466d31a37c029898.png
.. |image1| image:: /dirac/images/thumb/9/95/Dirac.png/41px-Dirac.png
.. |\\langle
i\|\\hat{A}\|i\\rangle| image:: /dirac/images/math/e/e/0/ee011c9a9fd1edd5b3dc6e3538f27934.png
.. |\\langle
i\|\\hat{A}\|j\\rangle| image:: /dirac/images/math/c/d/9/cd975f463a40ca67808a1655dd2715cb.png
.. |Powered by
MediaWiki| image:: /dirac/skins/common/images/poweredby_mediawiki_88x31.png
