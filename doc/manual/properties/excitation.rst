star(EXCITATION ENERGIES)

Calculate excitation energies using time dependent Hartree-Fock or DFT.
The excitation energies are found as the lowest generalized eigenvalues
of the electronic Hessian. DIRAC supports TDDFT kernels from all ground
state functionals included in the code. Currently the iterative
eigenvalue solver may fail to converge more than about twenty roots per
symmetry.

Define excitations and transition moments
=========================================

keyword(EXCITA)

::

    .EXCITA
    SYM N

Number of excitation energies N calculated in boson symmetry no. SYM.
This keyword can be repeated if you want excitation energies in more
than one boson symmetry.

keyword(OPERATOR)

Specification of a transition moment operator (see
:ref:`one_electron_operators` for details). This keyword can be given multiple
times to add more operators.

keyword(ANALYZE)

Analyze solution vectors and show the most important excitations at the
orbital level.

Control variational parameters
==============================

keyword(OCCUP)

For each fermion ircop give an :ref:`orbital_strings` of inactive orbitals from
which excitations are allowed. By default excitations from all occupied
orbitals are included in the generalized eigenvalue problem.

Example:

::

    .OCCUP
    1..3
    7,8

This would include excitations from gerade orbitals 1,2,3, and ungerade
orbitals 7 and 8.

keyword(VIRTUA)

For each fermion ircop give an :ref:`orbital_strings`
of virtual orbitals
to which excitations are allowed. By default excitations to all virtal
orbitals are included in the generalized eigenvalue problem.

keyword(SKIPEE)

Exclude all rotations between occupied positive-energy and virtual
positive-energy orbitals.

keyword(SKIPEP)

Exclude all rotations between occupied positive-energy and virtual
negative-energy orbitals.

Control reduced equations
=========================

keyword(MAXITR)

Maximum number of iterations.

*Default:*

::

    .MAXITR
     30

keyword(MAXRED)

Maximum dimension of matrix in reduced system.

*Default:*

::

    .MAXRED
     200

keyword(THRESH)

Threshold for convergence of reduced system.

*Default:*

::

    .THRESH
     1.0D-5

Control integral contributions
==============================

The user is encouraged to experiment with these options since they may
have an important effect on run time.

keyword(INTFLG)

Specify what two-electron integrals to include
(default: :ref:`.INTFLG` under :ref:`**HAMILTONIAN`).

keyword(CNVINT)

Set threshold for convergence before adding SL and SS integrals to
SCF-iterations.

*2 (real) Arguments:*

::

    .CNVINT
     CNVXQR(1) CNVXQR(2)

*Default:* Very large numbers.

keyword(ITRINT)

Set the number of iterations before adding SL and SS integrals to
SCF-iterations.

*Default:*

::

    .ITRINT
     1 1

Advanced/debug flags
====================

keyword(E2CHEK)

Generate a complete set of trial vector which implicitly allows the
explicit construction of the electronic Hessian. Only to be used for
small systems !

keyword(ONLYSF)

Only call FMOLI in sigmavector routine: only generate one-index
transformed Fock matrix  [Saue2003]_.

keyword(ONLYSG)

Only call FMOLI in sigmavector routine: 2-electron Fock matrices using
one-index transformed densities [Saue2003]_.

