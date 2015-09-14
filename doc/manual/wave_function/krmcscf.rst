:orphan:
 

star(KRMCSCF)

*2- and 4-component relativistic KR-MCSCF module*
=================================================
The original implementation written
by Joern Thyssen, Timo Fleig and Hans Joergen Aagaard Jensen,
parallelization and continuous improvement by Stefan Knecht and Hans
Joergen Aagaard Jensen. Linear symmetry adaptation by Stefan Knecht and Hans Joergen Aagaard Jensen.

For more details on the theory and actual implementation see :cite:`Jensen1996`, :cite:`Thyssen2004` , :cite:`Thyssen2008`, :cite:`Knecht2009`.


Mandatory keywords  under \*KRMCSCF
===================================

keyword(CI PROGRAM)

specifies the CI module behind KRMCSCF, choices are *LUCIAREL* and *GASCIP*. 
*Note*: GASCIP does not work with linear symmetry and for large-scale MCSCF (> :math:`1.0 \times 10^7` determinants) only
LUCIAREL can be used efficiently. *default*:

::

    .CI PROGRAM
    GASCIP

keyword(GAS SHELLS)

Specification of CI calculation with electron distribution in orbital (GAS) spaces. 
The first line contains the number of GA spaces to be used (1-7), 
followed by one line per GAS with a separation by a "/" of the min/max number of electrons in each GAS and the 
number of orbitals per fermion correp (either one (no inversion symmetry) or two (inversion
symmetry: *gerade* *ungerade*) entries per line). 

The first entry before the "/" gives the minimal number of accumulated (!) electrons
after consideration of this GAS, the second the corresponding maximum
number, separated by blanks. The minimum and maximum accumulated
occupations allow for a very flexible parameterization of the wave
function. All determinants fulfilling the occupation constraints will be
constructed. The second entry after the "/" gives  the 
number of orbitals per fermion correp. See also the open-shell input in :ref:`*SCF` which is similar to the syntax used
here. The design of a GAS scheme is non-trivial and should be motivated by the electronic structure of the system (e.g.
inner core, outer core, valence, virtual space). Sometimes it is useful to subdivide the valence space, for scientific reasons, or/and the
virtual space, for technical reasons (save core memory). See reference :cite:`Fleig2006a`,
pp. 27 for more details. 

Example for a CASSCF run with 10 electrons in 10 orbitals, also known as CAS(10,10):

::

    .GAS SHELLS
     1
     10 10 /  5  5

if all orbitals (fullMCSCF, feasible only for systems with 2-3 active electrons) should be included, use:

::

    .GAS SHELLS
     1
     2 2 / all


keyword(INACTIVE)

Inactive orbitals per fermion correp. *default*: all orbitals are active. The example below for a molecule with inversion
center marks the lowest 4 *gerade* and 2 *ungerade* orbitals as inactive, i.e. they are always kept doubly occupied in all determinants
of the CI expansion:

::

    .INACTIVE
     4 2

All remaining orbitals above the inactive + active space are considered as secondary orbitals, i.e. they are kept empty
in all determinants of the CI expansion. By default all orbital rotations between the inactive-active and active-secondary
space are included. 

Optional keywords under \*KRMCSCF
=================================

keyword(THRESH)

convergence threshold in the MCSCF gradient given as double precision
value. The present default value is quite tight and might be adapted if
one is only interested in MCSCF energies, e.g. :math:`1.d-03`. *default*:

::

    0.5d-05

keyword(SYMMETRY)

integer value giving the (boson/fermion) symmetry of the state to optimize on. *default*:

::

    1

If you run the calculation in **linear symmetry** you have to specify **2** :math:`\times \Omega` value
of the state to optimize on. The doubling stems from the fact
that we want to avoid non-integer input, e.g. in case of an odd number
of electrons we might have :math:`\Omega` = 1/2, 3/2, 5/2, etc. values and the
corresponding input for the :math:`\Omega` = 1/2 would then read as

::

    1

If we have a system with an even number of electrons and inversion
symmetry the input for an :math:`\Omega` = 2g state would read as

::

    4g

keyword(MAX MACRO)

integer value giving the maximum number of MACRO iterations in the MCSCF optimization. *default*:

::

    25

keyword(MAX MICRO)

integer value giving the maximum number of MICRO iterations in the MCSCF optimization. *default*:

::

    50

keyword(DELETE)

delete active-secondary e-e rotations (rotations between
electronic-electronic spinors) in the gradient and Hessian calculation
specified by an orbital string of virtual (secondary) orbitals for each
fermion correp. *default*: include all active-secondary e-e rotations. *example:* (deleting all e-e rotations for secondary orbitals 20,21 and
22 in fermion correp 1 (*gerade*) and 24,25,26,27 and 28 in fermion
correp 2 (*ungerade*)):

::

    .DELETE
    20,21,22
    24..28

keyword(SKIPEE)

skip e-e rotations (rotations between electronic-electronic spinors) in the gradient and Hessian calculation. *default*: include e-e rotations.

keyword(SKIPEP)

skip e-p rotations (rotations between electronic-positronic spinors) in the gradient and Hessian calculation. *default*: include e-p rotations.

keyword(PRINT)

raise default print level to the given integer value. Please use with
care as you may get millions of output lines if you choose a too high
value. *default*:

::

    .PRINT
     0

star(OPTIMI)

optional keywords under \*OPTIMI
================================

**NOTE: the following keywords must be placed under the input deck**

::

    *OPTIMI

keyword(NOOCCN)

compute natural orbital occupation numbers for the final optimized electronic state. *default*: do not compute natural orbital occupation numbers.


keyword(ANALYZ)

analyze the final CI wave function printing the coefficients for each determinant above a given threshold
:math:`10^{-2}`. *default*: do not analyze the final CI wave function.


keyword(MAX CI)

maximum number of initial CI iterations. *default*:

::

    .MAX CI
     5

keyword(MXCIVE)

maximum size of Davidson subspace. *default*: 3 times the number of eigenstates (see :ref:`KRCI_.CIROOTS`) to optimize on.
example:

::

    .MXCIVE
     3

