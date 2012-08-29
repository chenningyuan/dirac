HAMILTONIAN:  Sample inputs
===========================

Some examples of the Hamiltonian input.

Add dipole perturbation of strength 0.01 au in the z-direction to the
nonrelativistic Lévy-Leblond Hamiltonian:

::

    **HAMILTONIAN
    .LEVY-LEBLOND
    .OPERATOR
     ZDIPLEN
     COMFACTOR
     0.01

Do a 4-component B3LYP calculation with high precision in the radial
grid:

::

    **DIRAC
    .WAVE FUNCTIONS
    **WAVE FUNCTIONS
    .DHF
    **HAMILTONIAN
    .DFT
     B3LYP
    *DFT
    .RADINT
     1.0D-15

A functional can also be specified manually through the GGAKEY keyword.
In the next input example we use the PW86X and P86C functionals for
exchange and correlation respectively. We also tighten the angular grid
and make the pruning less pronounced:

::

    **DIRAC
    .WAVE FUNCTIONS
    **WAVE FUNCTIONS
    .DHF
    **HAMILTONIAN
    .DFT
     GGAKEY PW86X=1.0 P86C=1.0
    *DFT
    .ANGINT
     53
    .ANGMIN
     15

Setting various two-component Hamiltonians
-------------------------------------------

The user can select order of scalar and spin-orbit terms in the :ref:`.BSS` Hamiltonian.

Both scalar and spin-orbit relativistic effects up to infinite order:

::

    .BSS
     099

Scalar relativistic effects of type "from the beginning" up to the
infinite order, no spin-orbit interaction:

::

    .SPINFREE
    .BSS
     109

Scalar relativistic effects "from the end" up to the infinite order:

::

    .SPINFREE
    .BSS
     009

Traditional scalar relativistic "from the beginning" second-order
Douglas-Kroll-Hess Hamiltonian:

::

    .SPINFREE
    .BSS
     102

Scalar relativistic "from the end" second-order Douglas-Kroll-Hess Hamiltonian:

::

    .SPINFREE
    .BSS
     002

Douglas-Kroll-Hess Hamiltonian with first-order spin-orbit and second-order scalar relativistic effects:

::

    .BSS
     012

Douglas-Kroll-Hess Hamiltonian with second-order spin-orbit and second-order scalar relativistic effects:

::

    .BSS
     022


In the two-component variational scheme is possible to combine external
AMFI spin-orbit terms, [Ilias2001]_, together with BSS integrals - see the :ref:`\*AMFI` section. 
Note that AMFI provides only one-center atomic integrals.

Infinite order scalar and spin-orbit terms, (one-electron) mean-field
spin-same-orbit (MFSSO2) from AMFI:

::

    .BSS
    2999

This is in fact the 'maximum' of two-component relativity, resembling
Dirac-Coulomb Hamiltonian.

Second-order Douglas-Kroll-Hess spin-free from 'the beginning' with
first order spin-orbit (SO1) term plus mean-field spin-same-orbit
(MFSSO) from AMFI:

::

    .BSS
    2112

SO1 may come either from BSS-transformation or from AMFI. In the latter
case it is only for one-center and for point nucleus.

Infinite order scalar and spin-orbit terms, (one-electron) mean-field
spin-same and spin-other-orbit (MFSO2) from AMFI:

::

    .BSS
    3999

This mimics the Dirac-Coulomb-Gaunt Hamiltonian, as the spin-other-orbit
comes from the Gaunt interaction term.

DIRAC enables switching-off AMFI spin-orbit contributions from various
centers. See the keyword :ref:`.NOAMFC`.

Infinite order scalar terms "from the beginning", (one-electron) spin
orbit (SO1) and the mean-field spin-same-orbit (MFSSO2) terms
exclusively from AMFI:

::

    .BSS
    4109

Second order (DKH) scalar terms "from the beginning", (one-electron)
spin orbit (SO1) and the mean-field spin-same-orbit (MFSSO2) terms
exclusively from AMFI:

::

    .BSS
    4102

Infinite order scalar terms "from the end", (one-electron) spin orbit
(SO1) and the mean-field spin-same-orbit and spin-other-orbit (MFSO2)
terms from AMFI:

::

    .BSS
    5009

Infinite order scalar terms "from the end", and the (one-electron) spin
orbit term (SO1) from AMFI:

::

    .BSS
    6009

AMFI one-electron spin-orbit terms - SO1 - are currently for the point
nucleus.

For the BSS value 'axyz' of y=0, DIRAC employs spin-free picture change
transformation of property operators, although the system is not in the
boson (spin-free) symmetry for a>1.

Rough first order, DKH1 (not recommended for practical calculations):

::

    .BSS
    111

For comparison purposes only between BSS-SO1 and AMFI-SO1 atomic one-center integrals (point nucleus only) use

::

    .BSS
    001

