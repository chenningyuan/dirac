Manual:KR-CI
============

*2- and 4-component relativistic GASCI module* written by Timo Fleig
based on LUCIA by Jeppe Olsen, parallelization by Stefan Knecht

The KR-CI module is a string-based Hamiltonian-direct configuration
interaction (CI) program [Fleig2003]_, [Fleig2006]_ ,
based on the LUCIA code [Olsen1990]_.
It is capable of doing efficient
CI computations at arbitrary excitation level, e.g. FCI, SDCI, RASCI,
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
four-component formalism (see
`\*\*HAMILTONIAN </dirac/index.php/Manual:HAMILTONIAN>`_).

The technical limitations are roughly set by several 100 million
determinants in the CI expansion on PCs and common computing clusters
and several billions of determinants on supercomputers with ample
memory. For calculations involving more than 500.000 determinants and in
particular higher than double excitations it is strongly recommended to
use the KR-CI module in 2c- and 4c applications.

If desired, the program also computes 1-particle densities from
optimized CI wave functions from which the natural orbital occupations
are printed.

The program is also used for the computation of various one-electron
properties at the GASCI level [Knecht2009]_ .
Please refer to the `KR-CI property
module </dirac/index.php/Manual:KR-CI_properties>`_ page for more
information.

Mandatory keywords
==================

**NOTE: the following keywords must be placed under the input deck**

::

    *KRCICALC

.CI PROGRAM
-----------

LUCIAREL

*Default:* None.

.MK2REF
-------

DELTA\_MK2 reference value (integer value). Should be set to zero.
Change only if you know what you are doing.

*Default:* None.

.MK2DEL
-------

DELTA\_MK2 value (integer value). Should be set to the number of
explicitly correlated electrons in the calculation. Change only if you
know what you are doing.

*Default:* None.

.CIROOTS
--------

number of states in symmetry X to optimize on. This keyword may be
repeated several times for multi-root multi-symmetry calculations
step-by-step.

*Default:* None.

*example for 5 roots in (boson/fermion) irrep 1*

::

    .CIROOTS
    1 5

If you run the calculation in **linear symmetry** (**only available in
the development version**) you have to specify the **2 x OMEGA** value
of the state to optimize on. The doubling (*2 x*) stems from the fact
that we want to avoid non-integer input, e.g. in case of an odd number
of electrons we might have OMEGA=1/2, 3/2, 5/2, etc. values and the
corresponding input for one root with the OMEGA=1/2 would then read as

::

    .CIROOTS
    1 1

If we have a system with an even number of electrons and inversion
symmetry the input for two OMEGA=2g states read as

::

    .CIROOTS
    4g 2

.GASSH
------

Number and specification of GAS orbitals. Line with the number of GA
spaces used (1-7), followed by one line per GAS with number of orbitals
per fermion correp (either one (no inversion symmetry) or two (inversion
symmetry: g u) entries per line). The design of GAS is non-trivial and
should be motivated by the electronic structure of the system (e.g.
inner core, outer core, valence, virtual space). Sometimes it is useful
to subdivide the valence space, for scientific reasons, or/and the
virtual space, for technical reasons (save core memory). See reference [Fleig2006a]_,
pp. 27 for more details.

 *Default:* None.

.GASSPC
-------

specification of CI calculation. One line per GAS with 2 entries each:
The first entry gives the minimal number of accumulated (!) electrons
after consideration of this GAS, the second the corresponding maximum
number, separated by blanks. The minimum and maximum accumulated
occupations allow for a very flexible parameterization of the wave
function. All determinants fulfilling the occupation constraints will be
constructed.

*Default:* None.

optional keywords
=================

.MAX CI
-------

Number of CI iterations.

*Default:*

::

    .MAX CI
     5

.INACTIVE
---------

Inactive orbitals per fermion correp.

*Default:* All orbitals active.

.MXCIVE
-------

size of Davidson subspace.

example:

::

    .MXCIVE
     24

*Default:* 3 times the number of eigenstates to optimize on.

.RSTRCI
-------

*Default:* No restart. Optional restart from CI vector on file
KRCI\_CVECS.x where **x** is determined by the symmetry of the wave
function (symmetry == 1 --> **x = 1**; symmetry == 2 --> **x = 2**;
etc.). **NOTE:** in the DIRAC10 release version the **x** <==> symmetry
relation is: symmetry == 1 --> **x = a**; symmetry == 2 --> **x = b**;
etc.

::

    .RSTRCI
     1

*Default:* No restart.

.NOOCCN
-------

compute natural orbital occupation numbers for each electronic state.

*Default:*

::

    do not compute natural orbital occupation numbers.

.WRTFCK
-------

write Fock matrix to file KRMC\_FOCK where it can be re-read for a
possible restart of a calculation (see the following keyword).

*Default:*

::

    do not write Fock matrix to file.

.RDFOCK
-------

read Fock matrix from file KRMC\_FOCK.

*Default:*

::

    do not read Fock matrix from file.

.CHECKP
-------

enables a check point write of the current solution vectors to the file
KRCI\_CVECS.x (see above in .CIROOTS for an explanation of how **x** is
supposed to be replaced) during the Davidson iterations. A checkpoint
file will be written roughly every 6th iteration.

*Default:*

::

    do not write check points.

.ANALYZ
-------

analyze the final CI wave function printing the coefficients for each
determinant above a given threshold.

*Default:*

::

    do not analyze the final CI wave function.

.PRINT
------

raise default print level to the given integer value. Please use with
care as you may get millions of output lines if you choose a too high
value.

*Default:*

::

    print level is 0.

Parallel keywords
=================

.PARINT
-------

Specifies the availability of the two-electron integral file on each
co-workers scratch disk in parallel calculation.

*Default:* integral broadcast from master to co-workers

::

    .PARINT
     0

**Advanced options**

.IJKLRO
-------

enables the storage of the resorted two-electron integrals on file
IJKL\_REOD which can then be read-in in a restart. This avoids a
multiple reading from the 4IND\* files and subsequent resorting. The
keyword has to be present also in the restart to enable the potential
read-in procedure. Hint: This keyword may be combined with .MAX CI == 0
in a precedent step in order to save memory for the actual production
run. This is due to the fact that the resorting step itself requires the
in-core storage of the unsorted and sorted integrals.

*Default:*

::

    do not write the resorted integrals to file IJKL_REOD.

.IJKLSP
-------

enables the splitting of the resorted integrals among the co-workers
according to their needs in the computation of the sigma vector. It can
only be used in combination with.IJKLRO.

*Default:*

::

    do no split the resorted integrals among the co-workers.

Manual:KR-CI properties
=======================

*2- and 4-component relativistic KR-CI property module* written by
Stefan Knecht and Hans Joergen Aa. Jensen, parallelization by Stefan
Knecht

The KR-CI property module [Knecht2009]_
takes
advantage of the 2- and 4-component KR-CI module. It can be used for the
computation of:

::

    - permanent dipole moments in ground and excited states, 

::

    - transition dipole moments 

::

    - computation of ,  and  (Ω values) expectation values (only for linear molecules).

Upon request the analysis of other one-electron properties may also be
implemented. The module could in principle be used for each one-electron
operator that is specified in
`[1] <http://wiki.chem.vu.nl/dirac/index.php/Manual:One-electron_operators>`_.

**NOTE: the following keywords must be placed under the input deck**

::

    *KRCICALC

together with the appropriate CI calculation input (GASSH, GASSPC, etc.)

implemented keywords
====================

.DIPMOM
-------

compute the permanent dipole moments in electronic ground and excited
states.

.TRDM
-----

compute the transition dipole moments between electronic states.

.OMEGAQ
-------

compute the expectation values of the spin- and angular momentum
operator in z-direction and print the total
<\\hat{j}\_z> =
<\\hat{s}\_z> + <\\hat{l}\_z>| 
expectation value (**Ω** value) for each
electronic state.

::

    Note: this applies only to linear molecules.

