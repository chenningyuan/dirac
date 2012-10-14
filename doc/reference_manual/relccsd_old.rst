.. _relccsd:

=============================================
RELCCSD -- Coupled cluster module (old input)
=============================================

*Relativistic Coupled Cluster module RELCCSD* written by Lucas Visscher

This section specifies the input necessary for relativistic Coupled
Cluster calculations. The program is currently capable of doing energy
calculations at the MP2, CCSD, and CCSD(T) levels of theory [Visscher1996]_ , [Pernpointner2003]_  ,
while first order
properties are available only at the MP2 level of theory.

Open shell states can be calculated also via the Fock space Coupled
Cluster method [Visscher2001]_
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
non-iterative triple corrections [Visscher1996]_

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

Number of A operators in <<A; B>>.

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
