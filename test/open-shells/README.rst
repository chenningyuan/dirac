DIRAC open-shell tests
========================

Test dedicated to demonstrate Dirac HF-SCF open-shell capabilities.

The aim is to get converged ground and first excited states of 
several open-shell atoms and molecules.

For the sake of quick performance we use the STO-2G uncontracted basis sets and 
the 2-component (X2C+AMFI) Hamiltonian.

B 
--

The ground atomic state (^2P_1/2) is obtained by distributing 1 electron over 2 2p_12 spinors.
The reordering of valence spinors and distributing 1 electron over 4 2p_32 shells ensures converging to the first excited state ^2P_32.

C 
---

We obtained the ground state [He] 2s(2) 2p(2) by placing one open-shell electron into the 2p_12 shell and
the second into the 2p_32 spinor. The spinors SCF convergence with two open-shell levels 
and the proper degeneracy of the two 2p_32 shells was ensured by selecting proper ".OLEVEL" values.

N
---

The nitrogen atom...

O
--

For oxygen with four electrons in the 2p shell we show two open-shell configurations.
The first configuration of [He] 2s(2) 2p12(2) 2p32(2) has 2p12 shell doubly occupied (closed) 
and remaining two unpaired electrons are distributed over four  2p32 spinors.
The second configuration [He] 2s(2) 2p12(1) 2p32(3) has two open-shells with one unpaired electron averaged over
two 2p12 shells and three electron over four 2p32 shells.

 F
-----

The converged ground state X ^2P_3/2 is achieved with the help of the .OLEVEL keyword.
To get the first excited state ^2P_1/2 one resorts to the reordering of the valence molecular spinors,
where we obtain the desired occupation with swapped p12,p32 molecular spinors : 1s(2) 2s(2) 2p32(4) 2p12(1).

CH (^2Pi)
---------

The CH molecule...

