

.. _relccsd:

=======
RELCCSD
=======

The old namelist style input used up to DIRAC10 has been replaced by an input
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
multireference calculations. See further :ref:`CCFSPC`.

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

Number of active electrons. This variable determines the reference determinant to be used in the exponential expansion of the wave function. 
Since the default values correspond to the information passed on by the MOLTRA code on basis of the Hartree-Fock occupations 
and chosen range of active orbitals in MOLTRA, **there is usually no need to specify this variable manually**. 
Notable exception is the treatment of a high-spin open shell state with the single reference CCSD(T) ansatz. 
In such a case one first needs to determine the number of closed shell electrons in each irrep and 
then add to these occupations the open shell electrons (see below for an example). 

*Arguments:*

::

    Integer (NELEC(I),I=1,NFSYM*2).

*Default:*

::

     Active electrons in these irreps (written by **MOLTRA).

A simple example is the oxygen molecule in which the 1s electrons are kept frozen. The active space then consists of 12 electrons, 
of which there are 10 closed shell electrons. The closed shell electrons go in the nonbonding 2s :math:`\sigma_g` and :math:`\sigma_u` 
orbitals and the bonding 2p :math:`\sigma_g` and :math:`\pi_u` orbitals, while the 2 open shell electrons are distributed over the 
two antibonding 2p :math:`\pi_u` orbitals. Let us suppose that we want to take the :math:`M_S=1` state as our reference 
and assign alpha spin to both our open shell electrons. In the double group (:math:`D_{2h}^*`) we do not distinguish between 
:math:`\sigma` and :math:`\pi`, so we need to add all alpha electrons in gerade orbitals (2s :math:`\sigma_g` + 2p :math:`\sigma_g` 
+ 2p :math:`\pi_g` = 4), beta electrons in gerade orbitals (2s :math:`\sigma_g` + 2p :math:`\sigma_g` = 2), 
alpha electrons in ungerade orbitals (2s :math:`\sigma_u` + 2p :math:`\pi_u` = 3) 
and beta electrons in ungerade orbitals (2s :math:`\sigma_u` + 2p :math:`\pi_u` = 3). This gives the following occupation::

  .NELEC
  4 2 3 3

Note that this determinant does not represent the exact ground state of the oxygen molecule as this triplet is split by a few wave numbers 
due to spin-spin and (second order) spin-orbit interactions. The lowest state is the :math:`\Omega=0` component that cannot be 
represented by a single determinant. This state can be calculated using the Fock space coupled cluster method.

.NEL\_F1
--------

Number of electrons in the *gerade* irreps of the Abelian symmetry group. See below for more information.

.NEL\_F2
--------

Number of electrons in the *ungerade* irreps of the Abelian symmetry group. 

Coming back to the example given above, the oxygen molecule, we now show how this is done with the keywords given above. 
We again want to take the :math:`M_S=1` state as our reference. The irreps of :math:`D_{\infty h}^*` are ordered as 
1/2, -1/2, 3/2, -3/2,.... so we need to consider the :math:`\Omega` value (giving the projection on the molecular axis of both spin 
and orbital angular momentum) of the occupied oribitals. The :math:`\sigma` -orbitals go in the irreps 1/2 and -1/2 
while the :math:`\pi` -orbitals span the four irreps (1/2, -1/2, 3/2, -3/2). Putting an alpha electron in a :math:`\sigma` orbital will give 
an :math:`\Omega` value of 1/2, while putting it in a :math:`\pi` -orbital can either give -1/2 
(when put in the orbitals with orbital momentum -1) or 3/2 (when put in the +1 orbital). 
Similarly the beta electrons go in irreps -1/2 (for the :math:`\sigma`), -3/2 and 1/2 (for the :math:`\pi`). 
This makes the input for our example : ::

  .NEL_F1
  2 3 1 0
  .NEL_F2
  2 2 1 1 


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

low limit of orbital energies of active occupied
orbitals, which constitute the secondary Pi space. Could be
used in (1,0), (2,0) and (1,1) sectors. Arguments: real.


.EHMAX
~~~~~~

Maximum orbital energy of occupied orbitals forming the auxiliary (Pi)
space. Orbitals with energies higher than this energy are taken in the
primary (Pm) space.

This is upper limit of one-electronic energies of
active occupied orbitals, which constitute the secondary Pi space.
Could be used in (1,0), (2,0) and (1,1) sectors. Arguments: real.
 

.EPMIN
~~~~~~

Minimum orbital energy of virtual orbitals forming the auxiliary (Pi)
space. Orbitals with energies lower than this energy are taken in the
primary (Pm) space.

This is the low limit of orbital energies of active virtual
orbitals, which constitute the secondary Pi space.
Could be used in (0,1), (0,2) and (1,1) sectors. Arguments: real.
 

.EPMAX
~~~~~~

Maximum orbital energy of virtual orbitals forming the auxiliary (Pi)
space. Orbitals with energies higher than this energy are taken in the
secundary (Q) space and do not contribute to the model space.

This is the upper limit of one-electronic energies of active virtual
orbitals, which constitute the secondary Pi space. Could be
used in (0,1), (0,2) and (1,1) sectors. Arguments: real.
 

Other Intermediate Hamiltonian (IH) input parameters:
-----------------------------------------------------
 
For *experts* only.

Following keyowrds belong to the CCIH namelist section.
 
.IHSCHEME
~~~~~~~~~
Choose particular IH scheme. Arguments: Integer IHSCHEME =  1,  or 2.
 
The IHSCHEME =1 corresponds to the extrapolated IH
(XIH) approach, described in the paper [Eliav2005]_.

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
the paper [Landau2004]_.
 
Here the problematic denominators :math:`1/(E_{Pi}  -  E_{Q})` are
substituted simply by the factor 0.
 
Default: IHSCHEME = 2
 
Next key options are used only in case of XIH (IHSCHEME = 1).
 
.SHIFTH11  
~~~~~~~~~
Energy shift for the one-electronic excitations in (1,0) sector. Arguments: real.

.SHIFTH12 
~~~~~~~~~
Energy shift for the two-electronic excitations in (1,0) sector. Arguments: real.

.SHIFTH2     
~~~~~~~~
Energy shift for the two-electronic excitations in (2,0) sector. Arguments: real.

.SHIFTP11       
~~~~~~~~~
Energy shift for the one-electronic excitations in (0,1) sector. Arguments: real.

.SHIFTP12  
~~~~~~~~~
Energy shift for the two-electronic excitations in (0,1) sector. Arguments: real.

.SHIFTP2
~~~~~~~~
Energy shift for the two-electronic excitations in (0,2) sector. Arguments: real.
Usually we choose the approximate difference between the highest orbital energy belonging to Pi and the lowest
orbital energy belonging to the Pm space. 
Works only with the old style of RELCC input.

.SHIFTHP    
~~~~~~~~
Energy shift for the two-electronic excitations in (1,1) sector. Arguments: real

.AIH
~~~~
Compensation factor, used in expression (1). Arguments: real positive, not greater then 1.0.
 
.NIH 
~~~~
Compensation power, used in expression (1). Arguments: integer.
 
In the case of the limit: AIH=1.0 and NIH -> "infinity" ( NIH>100, in practice)
we have so called "full compensation" method, corresponding to the
extrapolation of the effective Hamiltonian from the intermediate one.
