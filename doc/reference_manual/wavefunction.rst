Manual:WAVE FUNCTION
====================

\*\*WAVE FUNCTION -- Get the wave function
==========================================

This section allows the specification of which wave function module(s)
to activate. By default no modules are activated. To activate any of
these modules you must also specify `.WAVE
FUNCTION </dirac/index.php/Manual:DIRAC#.WAVE_FUNCTION>`_ under
`\*\*DIRAC </dirac/index.php/Manual:DIRAC>`_, otherwise this input is
not read.

Note that the order below specifies the order in which the different
modules are called if you ask for more than one.

.SCF
----

Activates the Hartree-Fock/Kohn-Sham module.

Specification of the SCF module can be given in the
`\*SCF </dirac/index.php/Manual:WAVE_FUNCTION:SCF>`_ subsection.

If `.DFT </dirac/index.php/Manual:HAMILTONIAN#.DFT>`_ has been specified
under `\*\*HAMILTONIAN </dirac/index.php/Manual:HAMILTONIAN>`_, then a
Kohn-Sham calculation will be performed, otherwise a Hartree-Fock
calculation will be performed.

.RESOLVE
--------

Resolve open-shell states: do a small CI calculation to get the
individual energies of the states present in an
average-of-configurations open-shell Hartree-Fock calculation (see
`\*RESOLVE </dirac/index.php/Manual:WAVE_FUNCTION:RESOLVE>`_).

.MP2
----

Activates the `second-order Møller-Plesset
module </dirac/index.php/Manual:WAVE_FUNCTION:MP2CAL>`_.

.MVO
----

Calculate modified virtual orbitals (see
`\*MVOCAL </dirac/index.php/Manual:WAVE_FUNCTION:MVOCAL>`_).

.RELCCSD
--------

Activates the `coupled cluster
module </dirac/index.php/Manual:RELCCSD>`_ (and the
`\*\*MOLTRA </dirac/index.php/Manual:MOLTRA>`_ module to get 4-index
transformed integrals).

By default, molecular orbitals with orbital energy between -10 and +20
hartree (a.u.) are included, this can be modified in the
`\*\*MOLTRA </dirac/index.php/Manual:MOLTRA>`_ section.

.DIRRCI
-------

Activates the MOLFDIR CI module (and the
`\*\*MOLTRA </dirac/index.php/Manual:MOLTRA>`_ module to get 4-index
transformed integrals).

Specification of input for the MOLFDIR CI module is given in
`DIRRCI </dirac/index.php/Manual:DIRRCI>`_ and
`GOSCIP </dirac/index.php/Manual:GOSCIP>`_.

By default, molecular orbitals with orbital energy between -10 and +20
hartree (a.u.) are included, this can be modified in the
`\*\*MOLTRA </dirac/index.php/Manual:MOLTRA>`_ section.

.LUCITA
-------

Activates the `LUCITA CI module </dirac/index.php/Manual:LUCITA>`_ (and
the `\*\*MOLTRA </dirac/index.php/Manual:MOLTRA>`_ module to get 4-index
transformed integrals).

By default, molecular orbitals with orbital energy between -10 and +20
hartree (a.u.) are included, this can be modified in the
`\*\*MOLTRA </dirac/index.php/Manual:MOLTRA>`_ section.

**Pre-SCF orbital manipulations**

.REORDER MO'S
-------------

Interchange initial molecular orbitals prior to the SCF-calculation. The
start orbitals from DFCOEF are read and reordered.

For each fermion irrep give the new order of orbitals.

*Example:*

::

    .REORDER MO'S
    1..8,10,9

.ORBROT
-------

Jacobi rotations between pairs of orbitals.

On the line following the keyword, give first the rotation angle, then
on the following line(s) for each fermion irrep, give an `orbital
string </dirac/index.php/Manual:Orbital_strings>`_ of orbitals to
rotate.

**Post-SCF orbital manipulations**

.POST SCF REORDER MO'S
----------------------

Interchange converged molecular orbitals. The orbitals from DFCOEF are
read and reordered just before exiting the SCF subroutine.

For each fermion irrep give the new order of orbitals.

*Example:*

::

    .POST DHF REORDER MO'S
    1..8,10,9

.PHCOEF
-------

Phase adjustment of coefficients DFCOEF: make the largest element of a
given orbital real and positive.

.RELADC
-------

Activates the `propagator module </dirac/index.php/Manual:RELADC>`_ and
calculates the ionization spectrum by the (A)lgebraic (D)iagrammatic
(C)onstruction ADC. Also activates the
`\*\*MOLTRA </dirac/index.php/Manual:MOLTRA>`_ module to get 4-index
transformed integrals.

.KR CI
------

Activates the `2c/4c GAS-CI </dirac/index.php/Manual:LUCIAREL>`_ module
for the calculation of ground and excited states at the CI level.

.KRMCSCF
--------

Activates the `Kramers-restricted
MCSCF </dirac/index.php/Manual:KR-MCSCF>`_ module for the calculation of
ground and excited states at the CI level.

