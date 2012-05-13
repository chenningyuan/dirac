RELCCSD:Sample inputs
======================

Some examples of the namelist input used in DIRAC10.

Calculate the CCSD(T) energy and MP2 first order properties. The full
calculation is given as input example 1.energy+dipole in the test
directory:

::

     &RELCCSD TIMING=T, DOENER=T, DOFOPR=T &END
     &CCFOPR  DOMP2G=T &END


Oxygen molecule
-----------------
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

Coming back to the example given above, the oxygen molecule, we now
show how this is done with the keywords given in the reference manual. We again want to
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

:: 

 NELEC\_F1=2,3,1,0 NELEC\_F2=2,2,1,1


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

This example can be found as test #52 REFERENCE_IT!!! in the DIRAC test set.


Nitrogen atom
-------------

we give an
example for the high spin state of the nitrogen atom.... ADD ME !

