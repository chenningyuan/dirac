

star(KRCI) or star(KRCICALC) (deprecated)

*2- and 4-component relativistic GASCI module* written by Timo Fleig
based on LUCIA by Jeppe Olsen, parallelization by Stefan Knecht and Hans Joergen Aagaard Jensen.

The KRCI module is a string-based Hamiltonian-direct configuration
interaction (CI) program [Fleig2003]_, [Fleig2006]_ , based on the LUCIA code [Olsen1990]_.
It is capable of doing efficient CI computations at arbitrary excitation level, e.g. FCI, SDCI, RASCI,
and MRCI using general active spaces. The code is interfaced to
molecular integrals obtained in the relativistic 2c- and 4-component
framework and uses double-point group symmetry. It is implemented as a
full parallel version [Knecht2010a]_ .

A central feature of the program is the Generalized Active Space (GAS)
concept, in which the underlying total orbital space is subdivided into
a basically arbitrary number (save for an upper limit) of subspaces with
arbitrary occupation constraints. This is the most general approach to
orbital space subdivisions. The program uses DIRAC Kramers pairs from
either a closed- or an open-shell calculation in a relativistic two- and
four-component formalism (see :ref:`**HAMILTONIAN`).

The technical limitations are roughly set by several 100 million
determinants in the CI expansion on PCs and common computing clusters
and several billions of determinants on supercomputers with ample
memory.

If desired, the program also computes 1-particle densities from
optimized CI wave functions from which the natural orbital occupations
are printed.

The program can also be used for the computation of various one-electron
properties at the GASCI level [Knecht2009]_ .

**Mandatory keywords**
======================

keyword(CI PROGRAM)

specifies the CI module behind KRCI.

::

    .CI PROGRAM
    LUCIAREL

keyword(CIROOTS)

number of states (roots) in symmetry X to optimize on. This keyword may be
repeated several times for multi-root multi-symmetry calculations
step-by-step. *default:* none.

example for 5 roots in (boson/fermion) irrep 1:

::

    .CIROOTS
    1 5

If you run the calculation in **linear symmetry** you have to specify the **2 x OMEGA** value
of the state to optimize on. The doubling (*2 x*) stems from the fact
that we want to avoid non-integer input, e.g. in case of an odd number
of electrons we might have :math:`\Omega` = 1/2, 3/2, 5/2, etc. values and the
corresponding input for one root with the :math:`\Omega` =1/2 would then read as

::

    .CIROOTS
    1 1

If we have a system with an even number of electrons and inversion
symmetry the input for two :math:`\Omega` =2g states read as

::

    .CIROOTS
    4g 2

keyword(INACTI)

Inactive orbitals per boson symmetry, separated by commas. This keyword is not allowed in connection with the citype
choice :ref:`CITYPE` GASCI. *default*: all orbitals active. 
example for D:math:`_{2h}`:

::

    .INACTI
     1,1,0,0,1,0,0,2


keyword(GAS SHELLS)

::

    .GAS SHELLS
     4
      1  2 /  1  0
      6  8 /  1  2
     14 16 /  7  7
     16 16 / 30 24


**Optional keywords**
=====================

keyword(ANALYZ)

analyze the final CI wave function printing the coefficients for each determinant above a given threshold
:math:`10^{-2}`. *default*: do not analyze the final CI wave function.


keyword(MAXITR)

maximum number of CI iterations. *default*:

::

    .MAXITR
     100


keyword(RSTRCI)

Restart CI from vector(s) on file LUCVECT resp. LUCVECT.0 (parallel calculation):

::

    .RSTRCI
     1

*default*: no restart.

::

    .RSTRCI
     0

**Advanced options**
====================

keyword(IJKLRO)

enables the storage of the resorted two-electron integrals on file
IJKL\_REOD which can then be read-in in a restart. This avoids a
multiple reading from the 4IND\* files and subsequent resorting. The
keyword has to be present also in the restart to enable the potential
read-in procedure. Hint: This keyword may be combined with :ref:`MAX CI` == 0
in a precedent step in order to save memory for the actual production
run. This is due to the fact that the resorting step itself requires the
in-core storage of the unsorted and sorted integrals. *default*: do not write the resorted integrals to file IJKL_REOD.

keyword(IJKLSP)

enables the splitting of the resorted integrals among the co-workers
according to their needs in the computation of the sigma vector. It can
only be used in combination with :ref:`IJKLRO`. *default*: do no split the resorted integrals among the co-workers.

**KRCI properties**
===================

This part referes to the *2- and 4-component relativistic KR-CI property module* written by
Stefan Knecht and Hans Joergen Aa. Jensen, parallelization by Stefan Knecht

The KR-CI property module [Knecht2009]_ takes
advantage of the 2- and 4-component KRCI module. It can be used for the
computation of:

::

    - permanent dipole moments in ground and excited states,

::

    - transition dipole moments

::

    - computation of :math:`\hat{l}_z`, :math:`\hat{s}_z` and :math:`\hat{j}_z (\Omega)` expectation values (only for linear molecules).

Upon request the analysis of other one-electron properties may also be
implemented. The module could in principle be used for each one-electron
operator that is specified in :ref:`one_electron_operators`.

keyword(DIPMOM)

compute the permanent dipole moments in electronic ground and excited states.

keyword(TRDM)

compute the transition dipole moments between electronic states.

keyword(OMEGAQ)

*Note*: this keyword may only be used for linear molecules.

compute the expectation values of the spin- and angular momentum
operator in z-direction and print the total
:math:`\hat{j}_z` = :math:`\hat{s}_z` + :math:`\hat{l}_z` 
expectation value (:math:`\Omega` value) for each electronic state.

