
star(KRCI)

or (deprecated) star(KRCICALC)

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
number of orbitals per fermion correp. See also the open-shell input in :ref:`SCF` which is similar to the syntax used
here. The design of a GAS scheme is non-trivial and should be motivated by the electronic structure of the system (e.g.
inner core, outer core, valence, virtual space). Sometimes it is useful to subdivide the valence space, for scientific reasons, or/and the
virtual space, for technical reasons (save core memory). See reference [Fleig2006a]_,
pp. 27 for more details. example for 4 GAS spaces with an advanced excitation pattern:

::

    .GAS SHELLS
     4
      1  2 /  1  0
      6  8 /  1  2
     14 16 /  7  7
     16 16 / 30 24

if all remaining orbitals (fullCI) should be included in the last GAS space, use:

::

    .GAS SHELLS
     4
      1  2 /  1  0
      6  8 /  1  2
     14 16 /  7  7
     16 16 / all


keyword(INACTIVE)

Inactive orbitals per fermion correp. *default*: all orbitals active. example for 4 *gerade*  and 2 *ungerade* orbitals frozen in a molecule with inversion center:

::

    .INACTIVE
     4 2

**Optional keywords**
=====================

keyword(NOOCCN)

compute natural orbital occupation numbers for each electronic state. *default*: do not compute natural orbital occupation numbers.


keyword(ANALYZ)

analyze the final CI wave function printing the coefficients for each determinant above a given threshold
:math:`10^{-2}`. *default*: do not analyze the final CI wave function.


keyword(MAX CI)

maximum number of CI iterations. *default*:

::

    .MAX CI
     5

keyword(MXCIVE)

maximum size of Davidson subspace. *default*: 3 times the number of eigenstates (see :ref:`CIROOTS`) to optimize on.
example:

::

    .MXCIVE
     24

keyword(RSTRCI)

Restart CI from vector(s) on file KRCI\_CVECS.x where **x** is determined by the symmetry of the wave
function.

::

    .RSTRCI
     1

*default*: no restart.

::

    .RSTRCI
     0

*further infomation*: the convention in linear symmetry for KRCI_CVECS.x is the following where the offset/range for systems
with a gerade number of electrons is for x: offset --> 1; range: x =[1-64] and with an ungerade number of electrons x: offset -->
65; range: x =[65-128]. In more detail (with MJ == :math:`\Omega`):
    

a. for systems with a gerade number of e- and no inversion center:

::

    MJ =  0: KRCI_CVECS.1
    MJ = +1: KRCI_CVECS.2
    MJ = -1: KRCI_CVECS.3
    MJ = +2: KRCI_CVECS.4
    MJ = -2: KRCI_CVECS.5
    ...
    up to
    MJ = +32: KRCI_CVECS.64

    

b. for systems with a gerade number of e- and an inversion center:

::    

    MJ =  0g: KRCI_CVECS.1
    MJ = +1g: KRCI_CVECS.2
    MJ = -1g: KRCI_CVECS.3
    MJ = +2g: KRCI_CVECS.4
    MJ = -2g: KRCI_CVECS.5
    ...
    up to
    MJ = +16g: KRCI_CVECS.32

    and 
    MJ =  0u: KRCI_CVECS.33
    MJ = +1u: KRCI_CVECS.34
    MJ = -1u: KRCI_CVECS.35
    MJ = +2u: KRCI_CVECS.36
    MJ = -2u: KRCI_CVECS.37
    ...
    up to
    MJ = +16u: KRCI_CVECS.64

c. for systems with an ungerade number of e- and no inversion center:

::
    
    MJ = +1/2: KRCI_CVECS.65
    MJ = -1/2: KRCI_CVECS.66
    MJ = +3/2: KRCI_CVECS.67
    MJ = -3/2: KRCI_CVECS.68
    ...
    up to
    MJ = +63/2: KRCI_CVECS.128

d. for systems with an ungerade number of e- and an inversion center:

::    

    MJ = +1/2g: KRCI_CVECS.65
    MJ = -1/2g: KRCI_CVECS.66
    MJ = +3/2g: KRCI_CVECS.67
    MJ = -3/2g: KRCI_CVECS.68
    ...
    up to
    MJ = -31/2g: KRCI_CVECS.96

    and 

    MJ = +1/2u: KRCI_CVECS.97
    MJ = -1/2u: KRCI_CVECS.98
    MJ = +3/2u: KRCI_CVECS.99
    MJ = -3/2u: KRCI_CVECS.100
    ...
    up to
    MJ = -31/u: KRCI_CVECS.128


keyword(CHECKP)

enables a check point write of the current solution vectors to the file
KRCI\_CVECS.x (see above in :ref:`CIROOTS` for an explanation of how **x** is
supposed to be replaced) during the Davidson iterations. A checkpoint
file will be written roughly every 6th iteration. *default*: do not write check points.


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

