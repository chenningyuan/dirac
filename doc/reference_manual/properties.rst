PROPERTIES
===========

This section allows for the evaluation of a large number of molecular
properties. Available properties include:

-  Expectation values (e.g. dipole moment and electric field gradients).

-  Linear response properties (e.g. polarizability and NMR parameters).

-  Quadratic response properties (e.g. hyperpolarizabilities).

For convenience some common properties can be specified directly in this
section, which means that the user in principle does not need to know
how they are calculated. Note, however, that response functions are by
default static, but frequencies can be added in the relevant
subsections.

Properties which are not predefined must be specified in detail in the
relevant input section (see the `specification of one-electron
operators </dirac/index.php/Manual:One-electron_operators>`_).

By default no properties are calculated.

General control statements
==========================

.PRINT
------

Print level.

*Default:*

::

    .PRINT
     0

.ABUNDANCIES
------------

For properties that make reference to isotopes, give threshold level (in
% abundance) for isotopes to print.

*Default:*

::

    .ABUNDANCIES
     1.0

.RKBIMP
-------

Import coefficients calculated with restricted kinetic balance (RKB) in
a calculation using unrestricted kinetic balance (UKB). This option is a
simple way to generated restricted magnetic balance for the calculation
of NMR shieldings. This option works in the general SO case, but not in
the spinfree case since spinfree calculations are not possible with UKB.

Predefined electric properties
==============================

.DIPOLE
-------

Evaluate the electronic dipole moment (expectation value).

.QUADRUPOLE
-----------

Evaluate the electronic quadrupole moment (expectation value).

.EFG
----

Evaluate electric field gradients (expectation values).

Atomic centers may be restricted with
`.SELECT </dirac/index.php/Manual:INTEGRALS#.SELECT>`_ under
`\*\*INTEGRALS </dirac/index.php/Manual:INTEGRALS>`_.

.NQCC
-----

Evaluate nuclear quadrupole coupling constants (expectation values).

Atomic centers may be restricted with
`.SELECT </dirac/index.php/Manual:INTEGRALS#.SELECT>`_ under
`\*\*INTEGRALS </dirac/index.php/Manual:INTEGRALS>`_.

.POLARIZABILITY
---------------

Evaluate the static electronic dipole polarizability tensor (linear
response).

.FIRST ORDER HYPERPOLARIZABILITY
--------------------------------

Evaluate static electronic dipole first-order hyperpolarizability tensor
(quadratic response). Results are also given for the static electronic
dipole polarizability.

.VERDET
-------

Evaluate Verdet constants [Ekstrom2005]_
(quadratic response) for a dynamic electric field corresponding to Ruby
laser wavelength of 694 nm and a static magnetic field along the
propagation direction of the light beam (in this case, the default
frequencies of the quadratic response function thus become
ω\ :sub:`*B*`\  = 0.0656 and ω\ :sub:`*C*`\  = 0.0). A Verdet
calculation cannot be specified in combination with other quadratic
response calculations.

The frequencies can be changed using
`.BFREQ </dirac/index.php/Manual:PROPERTIES:QUADRATIC_RESPONSE#.BFREQ>`_
in `\*QUADRATIC
RESPONSE </dirac/index.php/Manual:PROPERTIES:QUADRATIC_RESPONSE>`_.

.TWO-PHOTON
-----------

Evaluate two-photon absorption cross sections (quadratic response). Give
the number of desired states in each boson symmetry. Cannot be specified
in combination with other quadratic response calculations.

*Example:* Point group with four boson irreps, (e.g.
*C*\ :sub:`2*v*`\ ).

::

    .TWO-PHOTON
     5 5 5 0

Predefined magnetic properties
==============================

.NMR
----

Evaluate nuclear magnetic shieldings and indirect spin-spin couplings
(linear response).

Atomic centers may be restricted with
`.SELECT </dirac/index.php/Manual:INTEGRALS#.SELECT>`_ under
`\*\*INTEGRALS </dirac/index.php/Manual:INTEGRALS>`_.

.SHIELDING
----------

Evaluate nuclear magnetic shieldings (linear response). `Print
level <#.PRINT>`_ 2 gives tensor and longer output. `Print
level <#.PRINT>`_ 4 gives the raw values in symmetry coordinates as
well.

Atomic centers may be restricted with
`.SELECT </dirac/index.php/Manual:INTEGRALS#.SELECT>`_ under
`\*\*INTEGRALS </dirac/index.php/Manual:INTEGRALS>`_.

.MAGNET
-------

Evaluate the magnetic susceptibilities tensor (linear response and
expectation values).

.SPIN-SPIN COUPLING
-------------------

Evaluate indirect spin-spin couplings (linear response).

Atomic centers may be restricted with
`.SELECT </dirac/index.php/Manual:INTEGRALS#.SELECT>`_ under
`\*\*INTEGRALS </dirac/index.php/Manual:INTEGRALS>`_.

.DSO
----

Evaluate the diamagnetic contribution to indirect spin-spin couplings as
an expectation value.

Atomic centers may be restricted with
`.SELECT </dirac/index.php/Manual:INTEGRALS#.SELECT>`_ under
`\*\*INTEGRALS </dirac/index.php/Manual:INTEGRALS>`_.

.NSTDIAMAGNETIC
---------------

Evaluate the diamagnetic contribution to nuclear magnetic shieldings as
an expectation value.

Atomic centers may be restricted with
`.SELECT </dirac/index.php/Manual:INTEGRALS#.SELECT>`_ under
`\*\*INTEGRALS </dirac/index.php/Manual:INTEGRALS>`_.

Other predefined properties
===========================

.MOLGRD
-------

Evaluate the molecular gradient, i.e. 

..|{\\frac{{\\partial E}}{{\\partial X\_A}}}|

where *X*\ :sub:`*A*`\  are the coordinates of the nuclei. This
is an expectation value of one- and two-electron operators. Normally the
molecular gradient evaluation is not invoked explicitly with this
keyword but rather implicitly in the geometry optimization module.

.PVC
----

Calculate matrix elements over the nuclear spin-independent
parity-violating operator, e.g. calculate energy differences between
enantiomers, see [Laerdahl1999]_ .

.RHONUC
-------

Calculate electronic density at the nucleus(nuclei).

Atomic centers may be restricted with
`.SELECT </dirac/index.php/Manual:INTEGRALS#.SELECT>`_ under
`\*\*INTEGRALS </dirac/index.php/Manual:INTEGRALS>`_.

.EFFDEN
-------

Calculate effective electronic density at the
nucleus(nuclei), [Knecht2011]_.

Atomic centers may be restricted with
`.SELECT </dirac/index.php/Manual:INTEGRALS#.SELECT>`_ under
`\*\*INTEGRALS </dirac/index.php/Manual:INTEGRALS>`_.


Manual:PROPERTIES:EXCITATION ENERGIES
=====================================

Calculate excitation energies using time dependent Hartree-Fock or DFT.
The excitation energies are found as the lowest generalized eigenvalues
of the electronic Hessian. DIRAC supports TDDFT kernels from all ground
state functionals included in the code. Currently the iterative
eigenvalue solver may fail to converge more than about twenty roots per
symmetry.

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
transformed Fock matrix  [Saue2003]_.

.ONLYSG
-------

Only call FMOLI in sigmavector routine: 2-electron Fock matrices using
one-index transformed densities [Saue2003]_.

\*EXPECTATION VALUE -- Evaluate expectation values
==================================================

This section gives directives for the calculation of expectation values.

.OPERATOR
---------

Specification of general one-electron operators.

See the `specification of one-electron
operators </dirac/index.php/Manual:One-electron_operators>`_ section for
more information and explicit examples.

**Advanced options**

.ORBANA
-------

Analyze contributions to an expectation value from the individual
orbitals.

.PROJEC
-------

Decompose an expectation value in contribution from fragment orbitals.
This action presupposes that a projection analysis has been performed,
see the `projection
analysis </dirac/index.php/Manual:ANALYZE:PROJECTION>`_ section for
details.

**Programmers options**

.PRINT
------

Print level.

*Default:*

::

    .PRINT
     0


Manual:PROPERTIES:LINEAR RESPONSE
=================================

Linear Response module written by T. Saue and H. J. Aa. Jensen [Saue2003]_.

The general form of a linear response function in the random phase
approximation (RPA) is [Saue2003]_, [Saue2002a]_

   :align: center
   :alt:  \\langle\\langle A ; B \\rangle\\rangle\_{\\omega} = -\\mathbf{E}\_A^{[1]\\dagger} \\bigl( E^{[2]} - \\omega S^{[2]} \\bigr)^{-1} \\mathbf{E}\_B^{[1]}
    \\langle\\langle A ; B \\rangle\\rangle\_{\\omega} = -\\mathbf{E}\_A^{[1]\\dagger} \\bigl( E^{[2]} - \\omega S^{[2]} \\bigr)^{-1} \\mathbf{E}\_B^{[1]}

where 
\|\\mathbf{E}\_A^{[1]}| 
and 
\|\\mathbf{E}\_B^{[1]}| 
are property
gradients corresponding to properties A and B, *E*\ :sup:`[2]`\  is the
molecular Hessian, *S*\ :sup:`[2]`\  is the metric, and ω is a
frequency.

The dimension of *E*\ :sup:`[2]`\  does in general not allow an explicit
inversion of the resolvent | \\bigl( E^{[2]} - \\omega S^{[2]}
\\bigr)^{-1} |

so one resorts to iterative techniques.

We write:

   :align: center
   :alt:  \\bigl( E^{[2]} - \\omega S^{[2]} \\bigr)^{-1} \\mathbf{E}\_B^{[1]} = \\mathbf{X}

    \\bigl( E^{[2]} - \\omega S^{[2]} \\bigr)^{-1} \\mathbf{E}\_B^{[1]}  = \\mathbf{X}

By a slight rearrangement we arrive at the linear response equations:

   :align: center
   :alt:  \\bigl( E^{[2]} - \\omega S^{[2]} \\bigr) \\mathbf{X} = \\mathbf{E}\_B^{[1]}

    \\bigl( E^{[2]} - \\omega S^{[2]} \\bigr) \\mathbf{X} = \\mathbf{E}\_B^{[1]}

The RPA equation is solved by expanding 
\|\\mathbf{X}| 
in a set of trial
vectors:

   :align: center
   :alt:  \\mathbf{X} = \\sum\_{i=1} \\mathbf{b}\_i a\_i; Y = [\\mathbf{b}\_1 \\mathbf{b}\_2 \\cdots \\mathbf{b}\_n]

    \\mathbf{X} = \\sum\_{i=1} \\mathbf{b}\_i a\_i; Y = [\\mathbf{b}\_1 \\mathbf{b}\_2 \\cdots \\mathbf{b}\_n]

and solving the reduced equation:

   :align: center
   :alt:  \\bigl( \\tilde{E}^{[2]} - \\omega \\tilde{S}^{[2]} \\bigr) \\mathbf{a} = \\tilde{\\mathbf{E}}\_B^{[1]}

    \\bigl( \\tilde{E}^{[2]} - \\omega \\tilde{S}^{[2]} \\bigr) \\mathbf{a} = \\tilde{\\mathbf{E}}\_B^{[1]}

in which appear: 
| \\tilde{E}^{[2]} = Y^\\dagger E^{[2]} Y;
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

   :align: center
   :alt: \\langle\\langle A ; B \\rangle\\rangle\_\\omega \\equiv \\langle\\langle B ; A \\rangle\\rangle\_\\omega

   \\langle\\langle A ; B \\rangle\\rangle\_\\omega \\equiv \\langle\\langle B ; A \\rangle\\rangle\_\\omega

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

Form all possible  \|\\langle\\langle A ; B \\rangle\\rangle\_\\omega|
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
efficiency of screening. CLARIFY!

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
transformed Fock matrix [Saue2003]_.

.ONLYSG
-------

Only call FMOLI in sigmavector routine: 2-electron Fock matrices using
one-index transformed densities [Saue2003]_.

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

\*MOLGRD -- Controls evaluation of the molecular gradient
=========================================================

The section gives directives for control of a single point evaluation of
the molecular gradient activated with
`.MOLGRD </dirac/index.php/Manual:PROPERTIES#.MOLGRD>`_ under
`\*\*PROPERTIES </dirac/index.php/Manual:PROPERTIES>`_.

The molecular gradient module is usually called from within the geometry
optimization module.


**Advanced options**

.TRICK
------

If activated the calculation of the SL and SS two-electron integral
contributions to the gradient is skipped if their contribution is
estimated to be small. An "empirical" estimate of the norms of the LS
and SS two-electron gradients based on the norm of the LL two-electron
gradient. If deactivated all contributions to gradient are calculated
and if activated then only "necessary" contributions are calculated.

*Default:* Deactivated.

.INTFLG
-------

Specify what two-electron integrals to include (see
`.INTFLG </dirac/index.php/Manual:HAMILTONIAN#.INTFLG>`_ under
`\*\*HAMILTONIAN </dirac/index.php/Manual:HAMILTONIAN>`_).

*Default:* `.INTFLG </dirac/index.php/Manual:HAMILTONIAN#.INTFLG>`_ from
`\*\*HAMILTONIAN </dirac/index.php/Manual:HAMILTONIAN>`_.

.SCREEN
-------

Screening threshold.

*Default:*

::

    .SCREEN
     1.0D-17

Which means that screening is turned off by default. This is due to an
enormous demand of memory.

**Programmers options**

.PRINT
------

Print level.

A print level of less than 3 gives only the total gradient. Print levels
from 3+ gives individual contributions to the gradient (kinetic energy
gradient, nuclear attraction gradient etc.). Print levels of 5+ print
some matrices, and 10+ give massive output!

*Default:* Print level from
`\*\*GENERAL </dirac/index.php/Manual:GENERAL>`_.

.NUMGRA
-------

Calculate the molecular gradient using numerical differentiation.

*Default:* Analytical evaluation if implemented.


\*NMR -- Controls evaluation of NMR parameters
==============================================

This section gives directives for the calculation of NMR parameters.

**Advanced options**

.LONDON
-------

Activate calculations of magnetic properties (NMR shielding constants
and magnetizabilities with London atomic orbitals.

*Default:* Use conventional atomic orbitals.

.GAUGEO
-------

Reads in a user defined gauge origin (in bohr). Note that both the gauge
origin and dipole origin are changed in order to have consistency.

*Default for conventional atomic orbitals:*

::

    .GAUGEO
     0.0 0.0 0.0

.GO\_ANG
--------

Same as `.GAUGEO </dirac/index.php/Manual:PROPERTIES:NMR#.GAUGEO>`_ but
reads coordinates in angstrom.

.USECM
------

Use the center of mass as the gauge origin.

.INTFLG
-------

Specify what two-electron integrals to include in the two-electron
London contributions to the magnetic field property gradient (see
`.INTFLG </dirac/index.php/Manual:HAMILTONIAN#.INTFLG>`_ under
`\*\*HAMILTONIAN </dirac/index.php/Manual:HAMILTONIAN>`_).

*Default:* `.INTFLG </dirac/index.php/Manual:HAMILTONIAN#.INTFLG>`_ from
`\*\*HAMILTONIAN </dirac/index.php/Manual:HAMILTONIAN>`_.

.NOTWO
------

Do not calculate the two-electron London contributions for the magnetic
field property gradient when London atomic orbitals are used.

.NOONEI
-------

Do not calculate the {H(0),T(B)} reorthonormalization terms for the
magnetic field property gradient when London atomic orbitals are used.

.NOORTH
-------

Do not calculate the {T(B),h(mK)} reorthonormalization contributions for
the expectation value term when London atomic orbitals are used.

.SYMCON
-------

Employ the symmetric connection for reorthonormalization terms when
using London atomic orbitals.

*Default:* Use the natural connection.


Manual:PROPERTIES:QUADRATIC RESPONSE
====================================

This section gives directives for the calculation of quadratic response
functions [Saue2002a]_, 

\|\\langle\\langle A ; B , C \\rangle\\rangle\_{\\omega\_1,\\omega\_2}|

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

Compute the transition moments \|\\langle i\|\\hat{A}\|i\\rangle|, 
where i goes from 1 to 5 in each symmetry (four symmetries in this case). The
zeros can be substituted for positive integers to generate elements \|\\langle i\|\\hat{A}\|j\\rangle|.

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

PROPERTIES:STEX
======================

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
files.\
 
please provide an
   example !

