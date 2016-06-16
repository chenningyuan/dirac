:orphan:
 

.. _relccsd:

============================================
Kramers unrestricted Coupled-Cluster methods
============================================

See :cite:`Visscher1996` for the initial CCSD(T) implementation. The Fock-space CC 
implementation is described in :cite:`Visscher2001`.

starstar(RELCC)

Specification of reference determinant, type of calculation, and general
settings.

keyword(FOCKSP)

Activate the Fock space module. This option should be used for
multireference calculations. See further :ref:`*CCFSPC`.
Because the first sector of Fock space gives the same result as a 
regular CCSD calculation, the latter calculation is switched off. If you
do wish to perform also a regular CC calculation (e.g. to get the CCSD(T) energy)
you need to activate this explicitly via :ref:`RELCC_.ENERGY` (see below).

keyword(ENERGY)

Activate the energy calculation. This is the preferred option for
calculations on closed shell or simple open shell systems and need not
be specified explictly in such cases. For Fock space calculations the
keyword switches on a separate single reference calculation done prior
to the FS calculation.

*Default:*

::

    Perform energy calculation.

keyword(GRADIENT)

Calculate the effective 1-particle density matrix. This option can be used to calculate
molecular properties. For closed-shell MP2 see :cite:`vanStralen2005`. For closed-shell CC 
see :cite:`Shee_CCgrad2016`.

*Default:*

::

    No gradient calculation.

keyword(NELEC)

Number of correlated electrons. This variable determines the reference determinant to be used in the exponential expansion of the wave function. 
Since the default values correspond to the information passed on by the MOLTRA code on basis of the Hartree-Fock occupations 
and chosen range of active orbitals in MOLTRA, **for CLOSED SHELL systems there is usually no need to specify this variable manually**. 
If you do chose to specify this manually, you should make sure to count electrons carefully as the numbers relate to the number of correlated
electrons rather than the total number. This number may therefore change if you change tresholds in the integral transformation.

For OPEN SHELL systems it is usually better to employ the multireference Fock space approach, except for simple cases such as 
a high-spin open shell state. In this case the single reference CCSD(T) ansatz works rater well. For such cases it is easier to use 
the :ref:`NELEC_OPEN` keyword that allows to just specify the distribution of open shell electrons over the irreps, but for backwards
compatibility we retain the older NELEC option as well.

*Arguments:*

::

    Integer (NELEC(I),I=1,NFSYM*2).

*Default:*

::

     Number of correlated electrons in closed shells in these irreps (written by **MOLTRA).

keyword(NELEC_OPEN)

Distribution of correlated open shell electrons over the irreps.  This determines the reference determinant for open shell single reference calculations.

This input should always be given for average-of-configuration SCF calculations that are followed by a CC calculation. The total number of electrons that is
given should correspond to the number of open shell electrons.

*Arguments:*

::

    Integer (NELEC_OPEN(I),I=1,NFSYM*2).

*Default:*

::

     Zero. (note that this default leads to wrong results because a closed shell ion will be calculated if no input is given).

keyword(NEL_F1)

Number of electrons in the *gerade* irreps of the Abelian symmetry group. This works like the NELEC keyword, but uses the 
supersymmetry possible for linear groups in which irreps are ordered as 1/2, -1/2, 3/2, -3/2, 5/2, ... (with the number
being the m_j value). 

keyword(NEL_F2)

Number of electrons in the *ungerade* irreps of the Abelian symmetry group. 

keyword(PRINT)

Print level.

*Default:*

::

    .PRINT
     0

keyword(TIMING)

Print detailed timing information.

*Default:*

::

    Only limited timing information is printed.

keyword(DEBUG)

Print debug information.

*Default:*

::

    Debug information is not printed.

star(CCENER)

Covers options related to energy.

keyword(NOMP2)

Deactivate MP2 calculation.

keyword(NOSD)

Deactivate CCSD calculation.

keyword(NOSDT)

Deactivate the calculation of perturbative triples. This is potentially
useful when running into memory problems for very big calculations and will also
save some CPU time.

keyword(MAXIT)

Set maximum number of iterations allowed to solve the CC equations.

keyword(MAXDIM)

Set maximum number of amplitude vectors used in the DIIS extrapolation.

keyword(NTOL)

Specify requested convergence (10^-NTOL) in the amplitudes.

keyword(NOSING)

Eliminate T1 amplitudes in the calculation (only interesting for test
purposes, this gives no computational speed-up).

keyword(NODOUB)

Eliminate T2 amplitudes in the calculation (only interesting for test
purposes, this gives no computational speed-up).
Deactivate contribution from doubles; corresponds to a CCS calculation.

star(CCFOPR)

Calculate first-order properties (expectation values) for the MP2 and CCSD wave
function.

keyword(MP2G)
Use MP2 wave function. 

keyword(CCSDG)
Use CCSD wave function.

keyword(NATORB)

Calculate natural orbitals (currently only for MP2 density matrix)

keyword(RELAXED)

Use orbital-relaxed density matrix (currently only for MP2).
The current default for MP2 and CC wave functions is to use the unrelaxed 
density matrix. This is computationally less expensive, but also less accurate.





star(CCFSPC)

Perform a Fock space MRCC calculation in which a model space is
correlated and then diagonalized to give CC energies for a set of
states.

keyword(DOIH)

Use the Intermediate Hamiltonian formalism in which an auxiliary space
is used to prevent the "intruder state" problem. Default: IH formalism
not used.

keyword(DOEA)

Calculate electron affinities (add one electron to the reference state,
allowing occupation of the active virtual orbitals)

keyword(DOIE)

Calculate ionization energies (remove one electron from the reference
state, allowing depletion of the active occupied orbitals)

keyword(DOEA2)

Calculate second electron affinities (add two electrons to the reference
state, allowing occupation of the active virtual orbitals)

keyword(DOIE2)

Calculate second ionization energies (remove two electrons from the
reference state, allowing depletion of the active occupied orbitals)

keyword(DOEXC)

Calculate excitation energies (allow excitation from the set of active
occupied orbitals to the set of active virtual orbitals)

keyword(NACTH)

Specification of the set of active hole orbitals (from which
ionization/excitation takes place)

keyword(NACTP)

Specification of the set of active particle orbitals (to which electron
attachment/excitation takes place)

keyword(MAXIT)

Maximum number of iterations allowed to solve the FSCC equations

keyword(MAXDIM)

Set maximum number of amplitude vectors used in the DIIS extrapolation.

keyword(NTOL)

Specify requested convergence (10^-NTOL) in the amplitudes.

keyword(GESTAT)

Specify the state number in the last active sector to pick the energy 
from (remember to account for degeneracies) for a state-specific FSCC 
geometry optimization based on a numerical gradient. 

star(CCIH)

Options for intermediate hamiltonian in FSCC.

keyword(EHMIN)

Minimum orbital energy of occupied orbitals forming the auxiliary (Pi)
space. Orbitals with energies lower than this energy are taken in the
secundary (Q) space and do not contribute to the model space.

low limit of orbital energies of active occupied
orbitals, which constitute the secondary Pi space. Could be
used in (1,0), (2,0) and (1,1) sectors. Arguments: real.

keyword(EHMAX)

Maximum orbital energy of occupied orbitals forming the auxiliary (Pi)
space. Orbitals with energies higher than this energy are taken in the
primary (Pm) space.

This is upper limit of one-electronic energies of
active occupied orbitals, which constitute the secondary Pi space.
Could be used in (1,0), (2,0) and (1,1) sectors. Arguments: real.

keyword(EPMIN)

Minimum orbital energy of virtual orbitals forming the auxiliary (Pi)
space. Orbitals with energies lower than this energy are taken in the
primary (Pm) space.

This is the low limit of orbital energies of active virtual
orbitals, which constitute the secondary Pi space.
Could be used in (0,1), (0,2) and (1,1) sectors. Arguments: real.

keyword(EPMAX)

Maximum orbital energy of virtual orbitals forming the auxiliary (Pi)
space. Orbitals with energies higher than this energy are taken in the
secundary (Q) space and do not contribute to the model space.

This is the upper limit of one-electronic energies of active virtual
orbitals, which constitute the secondary Pi space. Could be
used in (0,1), (0,2) and (1,1) sectors. Arguments: real.
 
Other Intermediate Hamiltonian (IH) input parameters
====================================================
 
For *experts* only.

Following keyowrds belong to the CCIH namelist section.
 
keyword(IHSCHEME)

Choose particular IH scheme. Arguments: Integer IHSCHEME =  1,  or 2.
 
The IHSCHEME =1 corresponds to the extrapolated IH
(XIH) approach, described in the paper :cite:`Eliav2005`.

Main idea: proper modification of the energetic denominators,
containing problematic Pi  -> Q transition.
The original denominator 1/(E_Pi  -  E_Q) , used during CC
iterations, is substituted by the following expression (1)

.. math::

 \frac{(1-[\frac{AIH*SHIFT}{(E_{Pi} - E_{Q} + SHIFT)}]^{NIH})}{\frac{(1-AIH*SHIFT}{(E_{Pi}  -  E_{Q} + SHIFT))}},

here AIH, SHIFT,NIH are parameters, specially chosen for
overcoming of the intruder states problem. These parameters
could be used in the procedure of the extrapolation of
the "exact" effective Hamiltonian solutions from
corresponding IH CC energies and wave functions.
 
IHSCHEME =2 corresponds to the simplified IH-2 approach, described in
the paper :cite:`Landau2004`.
 
Here the problematic denominators :math:`1/(E_{Pi}  -  E_{Q})` are
substituted simply by the factor 0.
 
Default: IHSCHEME = 2
 
Next key options are used only in case of XIH (IHSCHEME = 1).
 
keyword(SHIFTH11)

Energy shift for the one-electronic excitations in (1,0) sector. Arguments: real.

keyword(SHIFTH12)

Energy shift for the two-electronic excitations in (1,0) sector. Arguments: real.

keyword(SHIFTH2)

Energy shift for the two-electronic excitations in (2,0) sector. Arguments: real.

keyword(SHIFTP11)

Energy shift for the one-electronic excitations in (0,1) sector. Arguments: real.

keyword(SHIFTP12)

Energy shift for the two-electronic excitations in (0,1) sector. Arguments: real.

keyword(SHIFTP2)

Energy shift for the two-electronic excitations in (0,2) sector. Arguments: real.
Usually we choose the approximate difference between the highest orbital energy belonging to Pi and the lowest
orbital energy belonging to the Pm space. 
Works only with the old style of RELCC input.

keyword(SHIFTHP)

Energy shift for the two-electronic excitations in (1,1) sector. Arguments: real

keyword(AIH)

Compensation factor, used in expression (1). Arguments: real positive, not greater then 1.0.
 
keyword(NIH)

Compensation power, used in expression (1). Arguments: integer.
 
In the case of the limit: AIH=1.0 and NIH -> "infinity" ( NIH>100, in practice)
we have so called "full compensation" method, corresponding to the
extrapolation of the effective Hamiltonian from the intermediate one.

star(CCSORT)

Specialist options related to the sorting of two-electron integrals
and the calculation of the reference Fock matrix.

keyword(NORECMP)

Do not recompute the Fock matrix, but assume a diagonal matrix with the orbital energies 
taken from the SCF program on the dioagonal. This is usually not recommended as the latter
correspond to a restricted open shell expression and RELCCSD uses an unrestricted formalism.
For closed shell systems the two expressions are identical and this option merely suppresses
a build-in check on the accuracy of transformed integrals.

keyword(USEOE)

Ignore recomputed Fock matrix and use orbital energies supplied by the SCF program. This option is
sometimes useful for degenerate open shell cases in which case the perturbation theory for the unrestricted formalism
is not invariant for rotations among degenerate orbitals. It should only change the outcome
of the [T], (T) and -T energy corrections.
