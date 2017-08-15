:orphan:
 

.. _polprp:

===============================================================
Kramers restricted Polarization Propagator in the ADC framework
===============================================================

See :cite:`Pernpointner2014` for the initial POLPRP implementation. The transition
moment implementation and parallelization is described in :cite:`Pernpointner2017`.

starstar(POLPRP)

Specification of the overall calculation. Davidson diagonalizer (see DAVIDSON section)
is automatically invoked.

keyword(STATES)

Specifies the symmetries of the excited final states to be calculated. If you leave out this
keyword all symmetries are checked for the occurrence of possible final states and those are
calculated in the respective perturbational order. If the keyword is present, states have to be
specified according to the following input example

*Input example:*

::

   .STATES
   4
   1,3,17,19

This determines calculation for all eigenstates in the symmetries (irreducible representations) 1,3,17 and 19. 
Keep in mind that the numbering scheme for the double group ireps is a bit unconventional and depends on the
type of Hamiltonian used.


keyword(DOEXTE)

Presence of this keyword triggers an extended ADC(2) calculation in the sense that first-order contributions are
included in the 2p2h (satellite) block. This often helps to localize states having a sizeable amount
of double excitation character because those will alter their energy considerably when calculated in the extended scheme.




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
