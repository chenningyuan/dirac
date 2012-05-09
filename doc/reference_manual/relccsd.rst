Manual:RELCCSD
==============

*Relativistic Coupled Cluster module RELCCSD* written by Lucas Visscher

This section specifies the input necessary for relativistic Coupled
Cluster calculations. The program is currently capable of doing energy
calculations at the MP2, CCSD, and CCSD(T) levels of theory
[1]<#_note-visscher:relccsd1>`_`\ 
[2] <#_note-pernpointner:relccsd3>`_`\ , 

while first order
properties are available only at the MP2 level of theory.

Open shell states can be calculated also via the Fock space Coupled
Cluster method
[3] <#_note-visscher:relccsd2>`_`\  
in which
electrons are added to or subtracted from a closed shell core. All
sectors that involve no more than 2 replacements are available.

The input should be given in namelist form and is subdivided in a
general menu section and section describing energy, first-order,
second-order property and Fock-space Coupled Cluster calculations.

note of namelist: Some parts of DIRAC reads input in the Fortran *namelist* format. These
statements are put *after* the ``*END OF INPUT`` statement in the DIRAC
input file.

&RELCCSD -- Generic and control input
=====================================

Generic and control input.

DOENER
------

Calculate the CCSD(T) energy based on the standard single determinant
ansatz.

*Default:*

::

    DOENER=T

DOFSPC
------

Calculate CCSD energies using Fock space CC ansatz. In this method one
chooses a closed shell configuration from which the open shell
configurations of interest can be obtained by adding or subtracting one
or two electrons. The theory allows full coupling between the states
that belong to the same "sector" of Fock space (the set of
configurations that differ by the same amount of electrons from the
closed shell reference) and thus enables a multireference treatment of
open shell states. Energies are calculated and reported as electron
affinities (sectors [0,1] and [0,2]), ionization energies(sectors [1,0]
and [2,0]), or excitation energies (sector [1,1]) relative to the
reference closed shell state.

*Default:*

::

    DOFSPC=F

**Advanced options**

DOFOPR
------

Calculate the (relaxed) density matrix and first-order properties (as
specified in the section
`Manual:PROPERTIES </dirac/index.php/Manual:PROPERTIES>`_). These are
only available for the MP2 wave function. The MP2 density matrix is
**not** available for spinfree (or nonrelatvistic) calculations and the
program will stop when this combination of options is chosen.

*Default:*

::

    DOFOPR=F

NELEC
-----

Number of active electrons. This variable determines the reference
determinant to be used in the exponential expansion of the wave
function. Since the default values correspond to the information passed
on by the MOLTRA code on basis of the Hartree-Fock occupations and
chosen range of active orbitals in
`MOLTRA </dirac/index.php/Manual:MOLTRA#.ACTIVEMOLTRA>`_, there is
usually no need to specify this variable manually. Notable exception is
the treatment of a high-spin open shell state with the single reference
CCSD(T) ansatz. In such a case one first needs to determine the number
of closed shell electrons in each irrep and then add to these
occupations the open shell electrons (see below for an example).

*Arguments:* Integer (NELEC(I),I=1,NFSYM\*2).

*Default:* Active electrons in these irreps (written by
`\*\*MOLTRA </dirac/index.php/Manual:MOLTRA>`_).

 A simple example is the oxygen molecule in which the 1s electrons are
kept frozen. The active space then consists of 12 electrons, of which
there are 10 closed shell electrons. The closed shell electrons go in
the nonbonding 2s σ\ :sub:`g`\  and σ\ :sub:`u`\  orbitals and the
bonding 2p σ\ :sub:`g`\  and π\ :sub:`u`\  orbitals, while the 2 open
shell electrons are distributed over the two antibonding 2p
π\ :sub:`g`\  orbitals. Let us suppose that we want to take the
M\ :sub:`S`\ =1 state as our reference and assign alpha spin to both our
open shell electrons. In the double group (D\ :sub:`2h`\ \ :sup:`\*`\ )
we do not distinguish between σ and π, so we need to add all alpha
electrons in gerade orbitals (2s σ\ :sub:`g`\  + 2p σ\ :sub:`g`\  + 2p
π\ :sub:`g`\  = 4), beta electrons in gerade orbitals (2s σ\ :sub:`g`\ 
+ 2p σ\ :sub:`g`\  = 2), alpha electrons in ungerade orbitals (2s
σ\ :sub:`u`\  + 2p π\ :sub:`u`\ pi\_u = 3) and beta electrons in
ungerade orbitals (2s σ\ :sub:`u`\  + 2p π\ :sub:`u`\  = 3). This gives
the following occupation:

::

 NELEC=4,2,3,3

Note that this determinant does not represent the exact ground state of
the oxygen molecule as this triplet is split by a few wave numbers due
to spin-spin and (second order) spin-orbit interactions. The lowest
state is the Omega=0 component that cannot be represented by a single
determinant. This state can be calculated using the Fock space method
(see under `&CCFSPC </dirac/index.php/Manual:RELCCSD:CCFSPC#DOEA2>`_).

NELEC\_F1
---------

Number of electrons in the gerade irreps of the Abelian symmetry group.

*Arguments:* Integer (NELEC\_F1(I),I=1,number of gerade irreps).

*Default:* Values set by `\*\*MOLTRA </dirac/index.php/Manual:MOLTRA>`_
or `NELEC <#NELEC>`_.

NELEC\_F2
---------

Number of electrons in the ungerade irreps of the Abelian symmetry
group.

*Arguments:* Integer (NELEC\_F2(I),I=1,number of ungerade irreps).

*Default:* Values set by `\*\*MOLTRA </dirac/index.php/Manual:MOLTRA>`_
or `NELEC <#NELEC>`_.

 Coming back to the example given above, the oxygen molecule, we now
show how this is done with the keywords given above. We again want to
take the M\_S=1 state as our reference. The irreps of Dinfh\* are
ordered as 1/2,-1/2,3/2,-3/2,.... so we need to consider the omega value
(giving the projection on the molecular axis of both spin and orbital
angular momentum) of the occupied oribitals. The sigma orbitals go in
the irreps 1/2 and -1/2 while the pi-orbitals span the four irreps
(1/2,-1/2,3/2,-3/2). Putting an alpha electron in a sigma orbital will
give an omega value of 1/2, while putting it in a pi-orbital can either
give -1/2 (when put in the orbitals with orbital momentum -1) or 3/2
(when put in the +1 orbital). Similarly the beta electrons go in irreps
-1/2 (for the sigma), -3/2 and 1/2 (for the pi). This makes the input
for our example :

 NELEC\_F1=2,3,1,0 NELEC\_F2=2,2,1,1

IPRNT
-----

Print level. Level 1 gives information on the convergence of the CCSD
iterations and may be useful.

*Default:*

::

    IPRNT=0

TIMING
------

Write timing information.

*Default:*

::

    TIMING=F

The detailed `NELEC\_F1 <#NELEC_F1>`_ and `NELEC\_F2 <#NELEC_F2>`_ input
is mainly of interest in atomic calculations where the program may use
the D\_{\\infty h}| symmetry group. In that case this option permits
precise specification of the reference determinant.
`Here </dirac/index.php/Manual:RELCCSD:Sample_inputs>`_ we give an
example for the high spin state of the nitrogen atom.

**Programmers options**

DEBUG
-----

Write debug information.

*Default:*

::

    DEBUG=F

DOSOPR
------

Calculate second-order properties (only for testing and based on the HF
wave function).

*Default:*

::

    DOSOPR=F

DOSORT
------

Sort the incoming integrals over molecular spinors. This is a necessary
step that can only be skipped in special (restart) cases.

*Default:*

::

    DOSORT=T

NFROZ
-----

Number of frozen electrons in each irrep. Makes T1 and T2 amplitudes
that concern frozen spinors zero. Only for test purposes.

*Arguments:* Integer (NFROZ(I),I=1,NFSYM\*2).

*Default:* 0

INTERFACE
---------

This keyword controls how RELCCSD should handle the (transformed)
integral files coming from `MOLTRA </dirac/index.php/Manual:MOLTRA>`_.
Used together with scheme 6 this keyword should be specified as

::

    INTERFACE='DIRAC6'

Since scheme 6 is the default for MOLTRA and RELCCSD, this is normally
not necessary. However, if one uses scheme 4 in MOLTRA, the keyword
**MUST** be specified as shown below.

::

    INTERFACE='DIRAC'

&CCENER -- Energy module
========================

Energy module.

DOMP2
-----

Calculate the MP2 energy.

*Default:*

::

    DOMP2=T

DOCCSD
------

Calculate the CCSD energy.

*Default:*

::

    DOCCSD=T

DOCCSDT
-------

Calculate the CCSD(T) energy. The program calculates three different
non-iterative triple corrections
.. \ :sup:``[1] <#_note-visscher:relccsd1>`_`\ .

*Default:*

::

    DOCCSDT=T

**Advanced options**

NTOL
----

Convergence tolerance of the DIIS error vector (10\ :sup:`--NTOL`\ ).

*Default:*

::

    NTOL=12

MAXIT
-----

Maximum number of DIIS iterations to solve the CCSD equations.

*Default:*

::

    MAXIT=30

MAXDIM
------

Maximum dimension of the reduced space used in the DIIS extrapolation.

*Default:*

::

    MAXDIM=8

**Programmers options**

NOCCS
-----

Ignore single excitations so that the CCD formalism is used. Only for
test purposes.

*Default:*

::

    NOCCS=F

NOCCD
-----

Ignore double excitations so that the CCS formalism is used. Only for
test purposes.

*Default:*

::

    NOCCD=F

&CCSORT -- Sorting module
=========================

Sorting module.

USEOE
-----

Use orbital energies from the HF calculation instead of the diagonal
elements of the reconstructed Fock matrix. Putting USEOE=F may solve
problems with low-symmetry open shell systems, in which a Kamers pair is
split over the occupied and virtual space, giving an orbital energy
difference of zero in the default definition.

*Default:*

::

    USEOE=T

NORECMP
-------

Do not recompute the Fock matrix for the determinant selected as
reference for the coupled cluster expansion. This is usually done as a
check on the numerical accuracy of the MO transformation. If you specify
both USEOE and NORECMP the reference Hamiltonian is given by the
original orbital energies (taken from the SCF program) and the
transformed 2-electron integrals. This makes it possible to apply the
molecular mean field approximation.

*Default:*

::

    NORECMP=F

For the open shell calculation, turn off this flag. Use it with USEOE
option.

::

    &RELCCSD NELEC=4,3 TIMING=T, IPRNT=1, DOSORT=T, DOENER=T &END
    &CCSORT USEOE=F, NORECMP=F &END
    &CCENER DOMP2=T, DOCCSD=T, DOCCSDT=F &END



&CCFOPR -- First-order properties
=================================

First-order properties.

Only MP2 first order properties are currently available. They are
computed by forming a relaxed density matrix, which is done (partly) in
AO basis to obtain all relaxation contributions. The properties that are
evaluated are those specified in the expectation value part of the
input.

DOMP2G
------

Calculate the MP2 gradient, that is, expectation values using the
relaxed one-electron density matrix. Note that this is not the molecular
gradient, but other first order properties as specified in the
expectation value part of the input.

*Default:*

::

    DOMP2G=F

**Advanced options for the old formalism**

DOMP2GOLD
---------

Calculate the MP2 gradient via the old (approximate) MO formalism.

*Default:*

::

    DOMP2GOLD=F

NEOPER
------

Number of first order property operators \\langle E \\rangle| (only
relevant for the old formalism).

*Default:*

::

    NEOPER=1

NAMEE
-----

Names of *E* operators in image2|.

*Arguments:* Character (A8) (NAMEE(I),I = 1,NEOPER).

*Default:*

::

    NAMEE='ZDIPLEN'

NTOL
----

Convergence tolerance of the DIIS error vector (10\ :sup:`--NTOL`\ ).

*Default:*

::

    NTOL=12

MAXIT
-----

Maximum number of iterations to solve the Z-vector equations.

*Default:*

::

    MAXIT=30

MAXDIM
------

Maximum dimension of the reduced space.

*Default:*

::

    MAXDIM=8


&CCSOPR -- Second-order properties
==================================

Second-order properties

Only the SCF wave function can be used at the moment (RPA). The names of
the property integrals should first be defined in the
`.PRPTRA </dirac/index.php/Manual:MOLTRA#.PRPTRA>`_ section of
`\*\*MOLTRA </dirac/index.php/Manual:MOLTRA>`_.

NAOPER
------

Number of A operators in \\langle\\langle A ; B
\\rangle\\rangle\_\\omega|.

*Default:*

::

    NAOPER=1

NBOPER
------

Number of B operators in image1|.

*Default:*

::

    NBOPER=1

NFREQ
-----

Number of frequencies in image2|.

*Default:*

::

    NFREQ=1

NAMEA
-----

Names of A operators in image3|.

*Arguments:* Character (A8) (NAMEA(I),I = 1,NAOPER).

*Default:*

::

    NAMEA='ZDIPLEN'

NAMEB
-----

Names of B operators in image4|.

*Arguments:* Character (A8) (NAMEB(I),I = 1,NBOPER).

*Default:*

::

    NAMEB='ZDIPLEN'

EFREQ
-----

Values of frequencies in image5|.

*Arguments:* Real (EFREQ(I),I = 1,NFREQ).

*Default:*

::

    EFREQ=0.0

**Advanced options**

NTOL
----

Convergence tolerance of the DIIS error vector (10\ :sup:`--NTOL`\ ).

*Default:*

::

    NTOL=12

MAXIT
-----

Maximum number of iterations to solve the Z-vector equations.

*Default:*

::

    MAXIT=30

MAXDIM
------

Maximum dimension of the reduced space.

*Default:*

::

    MAXDIM=8


&CCFSPC -- Fock space module
============================

Fock space Coupled Cluster module.

DOEXC
-----

Calculate excitation energies.

*Default:*

::

    DOEXC=F

DOEA
----

Calculate electron affinities.

*Default:*

::

    DOEA=F

DOIE
----

Calculate ionization energies.

*Default:*

::

    DOIE=F

DOEA2
-----

Calculate electron affinities (add two electrons).

*Default:*

::

    DOEA2=F

DOIE2
-----

Calculate ionization energies (remove two electrons).

*Default:*

::

    DOIE2=F

NACTP
-----

Number of active spinors in which electrons may be put. This variable
determines the active space for the electron affinities (EA)
calculation.

*Arguments:* Integer (NACTP(I),I=1,NFSYM\*2).

*Default:* No default; should always be set when using `DOEA <#DOEA>`_
or `DOEA2 <#DOEA2>`_!

NACTH
-----

Number of active spinors in which holes can be created. This variable
determines the active space for the ionization energies (IE)
calculation.

*Arguments:* Integer (NACTH(I),I=1,NFSYM\*2).

*Default:* No default; should always be set when using `DOIE <#DOIE>`_
or `DOIE2 <#DOIE2>`_!

Example
---------

An example of the use of the Fock space method concerns the calculation
of the lowest three states of molecular oxygen. These are obtained by
distributing two electrons in the degenerate π\ :sub:`g`\  orbitals. The
X ground state is a triplet that is split by second-order spin-orbit
coupling and spin-spin interactions into a lowest M\ :sub:`S`\ =0
component and two higher M\ :sub:`S`\ =±2 components. This case is also
discussed under the RELCCSD input as an example of the use of the
`NELEC </dirac/index.php/Manual:RELCCSD:RELCCSD#NELEC>`_ keyword. In the
Fock space approach we first perform a calculation on
O\ :sub:`2`\ \ :sup:`+2`\  and calculate the three X states, the
degenerate open shell singlet a and the singlet b state in one single
Fock space run. The FSCCSD input reads then

:: 

 &RELCCSD NELEC=2,2,3,3, DOFSPC=T, DOENER=F &END

 &CCFSPC MAXIT=100, NACTP=2,2,0,0, DOEA2=T, &END

This example can be found as test #52 REFERENCE_IT in the Dirac test set.

**Advanced options**

NTOL
----

Convergence tolerance of the DIIS error vector (10\ :sup:`--NTOL`\ ).

*Default:*

::

    NTOL=12

MAXIT
-----

Maximum number of iterations to solve the CCSD equations.

*Default:*

::

    MAXIT=30

MAXITmn
~~~~~~~

Maximum number of iterations to solve the CCSD equations on each of the
Fock-space sectors.

*Default:*

::

    MAXIT00 = MAXIT01 = MAXIT10 = MAXIT11 = MAXIT20 = MAXIT02 = MAXIT

MAXDIM
------

Maximum dimension of the reduced space.

*Default:*

::

    MAXDIM=8

**Programmers options**

FSSECT
------

Active Fock space sectors (m,n) specifying the number of valence holes
(m) and valence particles (n) with respect to the reference determinant
(ordered 00,01,10,11,02,20).

*Default:*

::

    FSSECT=0,0,0,0,0,0


&CCIH -- Intermediate Hamiltonian Fock space CC functionality
=============================================================

EHMIN
-----

Lower boundary for the part of the *P*\ :sub:`*I*`\  space associated
with orbitals in which holes can be created (occupied orbitals for
sector 00).

*Default:*

::

    EHMIN = 1000.

EHMAX
-----

Upper boundary for the part of the *P*\ :sub:`*I*`\  space associated
with orbitals in which holes can be created (occupied orbitals for
sector 00).

*Default:*

::

    EHMAX = 100.

EPMIN
-----

Lower boundary for the part of the *P*\ :sub:`*I*`\  space associated
with orbitals in which particles can be created (virtual orbitals for
sector 00).

*Default:*

::

    EPMIN = -1.

EPMAX
-----

Upper boundary for the part of the *P*\ :sub:`*I*`\  space associated
with orbitals in which particles can be created (virtual orbitals for
sector 00).

*Default:*

::

    EPMAX = -1000.


Manual:RELCCSD:Sample inputs
============================

Some examples of the namelist input used in DIRAC10.

Calculate the CCSD(T) energy and MP2 first order properties. The full
calculation is given as input example 1.energy+dipole in the test
directory:

::

     &RELCCSD TIMING=T, DOENER=T, DOFOPR=T &END
     &CCFOPR  DOMP2G=T &END

Manual:RELCCSD:New input
========================

The old namelist style input (described above) used in DIRAC10 will be
replaced by an input that is more consistent with the rest of the input
in the next release. This input is briefly described for those of you
who have access to the developer version of DIRAC.

The old namelist style input used up to DIRAC10 is replaced by an input
that is more consistent with the rest of the input. In DIRAC11 it is
possible to use both styles of inputs but we have the intention to phase
out the namelist input in upcoming releases.

\*\*RELCC
=========

Specification of reference determinant, type of calculation, and general
settings.

.FOCKSP
-------

Activate the Fock space module. This option should be used for
multireference calculations. See further
`CCFSPC </dirac/index.php/Manual:RELCCSD:CCFSPC>`_

.ENERGY
-------

Activate the energy calculation. This is the preferred option for
calculations on closed shell or simple open shell systems.

*Default:*

::

    No energy calculation.

.GRADIENT
---------

Calculate the effective 1-particle density matrix (currently only for a
closed shell MP2 wave function). This option can be used to calculate
molecular properties.

*Default:*

::

    No gradient calculation.

.NELEC
------

Specify the occupation of the reference determinant (both for CC and
FSCC calculations). More information can be found here:
`NELEC </dirac/index.php/Manual:RELCCSD:RELCCSD>`_

.NEL\_F1
--------

Number of electrons in the gerade irreps of the Abelian symmetry group.
See also `NELEC\_F1 </dirac/index.php/Manual:RELCCSD:RELCCSD>`_

.NEL\_F2
--------

Number of electrons in the ungerade irreps of the Abelian symmetry
group. See also `NELEC\_F2 </dirac/index.php/Manual:RELCCSD:RELCCSD>`_

.PRINT
------

Print level.

*Default:*

::

    .PRINT
     0

.DEBUG
------

Print debug information.

*Default:*

::

    Debug information is not printed.

\*CCENER
--------

Covers options related to energy.

.DOMP2
~~~~~~

Calculate the MP2 energy.

.MAXIT
~~~~~~

Set maximum number of iterations allowed to solve the CC equations.

.MAXDIM
~~~~~~~

Set maximum number of amplitude vectors used in the DIIS extrapolation.

.NTOL
~~~~~

Specify requested convergence (10^-NTOL) in the amplitudes.

.NOSING
~~~~~~~

Eliminate T1 amplitudes in the calculation (only interesting for test
purposes, this gives no computational speed-up).

.NODOUB
~~~~~~~

Eliminate T2 amplitudes in the calculation (only interesting for test
purposes, this gives no computational speed-up).

\*CCFOPR
--------

Calculate first-order properties (expectation values) for the MP2 wave
function.

\*CCFSPC
--------

Perform a Fock space MRCC calculation in which a model space is
correlated and then diagonalized to give CC energies for a set of
states.

.DOIH
~~~~~

Use the Intermediate Hamiltonian formalism in which an auxiliary space
is used to prevent the "intruder state" problem. Default: IH formalism
not used.

.DOEA
~~~~~

Calculate electron affinities (add one electron to the reference state,
allowing occupation of the active virtual orbitals)

.DOIE
~~~~~

Calculate ionization energies (remove one electron from the reference
state, allowing depletion of the active occupied orbitals)

.DOEA2
~~~~~~

Calculate second electron affinities (add two electrons to the reference
state, allowing occupation of the active virtual orbitals)

.DOIE2
~~~~~~

Calculate second ionization energies (remove two electrons from the
reference state, allowing depletion of the active occupied orbitals)

.DOEXC
~~~~~~

Calculate excitation energies (allow excitation from the set of active
occupied orbitals to the set of active virtual orbitals)

.NACTH
~~~~~~

Specification of the set of active hole orbitals (from which
ionization/excitation takes place)

.NACTP
~~~~~~

Specification of the set of active particle orbitals (to which electron
attachment/excitation takes place)

.MAXIT
~~~~~~

Maximum number of iterations allowed to solve the FSCC equations

.MAXDIM
~~~~~~~

Set maximum number of amplitude vectors used in the DIIS extrapolation.

.NTOL
~~~~~

Specify requested convergence (10^-NTOL) in the amplitudes.

\*CCIH
------

Options for intermediate hamiltonian in FSCC.

.EHMIN
~~~~~~

Minimum orbital energy of occupied orbitals forming the auxiliary (Pi)
space. Orbitals with energies lower than this energy are taken in the
secundary (Q) space and do not contribute to the model space.

.EHMAX
~~~~~~

Maximum orbital energy of occupied orbitals forming the auxiliary (Pi)
space. Orbitals with energies higher than this energy are taken in the
primary (Pm) space.

.EPMIN
~~~~~~

Minimum orbital energy of virtual orbitals forming the auxiliary (Pi)
space. Orbitals with energies lower than this energy are taken in the
primary (Pm) space.

.EPMAX
~~~~~~

Maximum orbital energy of virtual orbitals forming the auxiliary (Pi)
space. Orbitals with energies higher than this energy are taken in the
secundary (Q) space and do not contribute to the model space.



References
----------

#. `↑ <#_ref-visscher:relccsd1_0>`_ L. Visscher, T. J. Lee, and K. G.
   Dyall, *Formulation and implementation of a relativistic unrestricted
   coupled cluster method including noniterative connected triples*, J.
   Chem. Phys. **105**, *8769* (1996) `electronic
   version <http://dx.doi.org/10.1063/1.472655>`_.
#. `↑ <#_ref-pernpointner:relccsd3_0>`_ M. Pernpointner and L. Visscher,
   *Parallelization of four-component calculations. II. Symmetry- driven
   parallelization of the 4-spinor CCSD algorithm*, J. Comput. Chem.
   **24**, *754* (2003) `electronic
   version <http://dx.doi.org/10.1002/jcc.10215>`_.
#. `↑ <#_ref-visscher:relccsd2_0>`_ L. Visscher, E. Eliav, and U.
   Kaldor, *Formulation and implementation of the relativistic
   Fock-space coupled cluster method for molecules*, J. Chem. Phys.
   **115**, *9720* (2001) `electronic
   version <http://dx.doi.org/10.1063/1.1415746>`_.

