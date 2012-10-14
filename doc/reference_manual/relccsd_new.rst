

=============================================
RELCCSD -- Coupled cluster module (new input)
=============================================

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

