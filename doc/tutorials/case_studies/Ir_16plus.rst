==============================
Getting excited states Ir(16+)
==============================

In the following we are interested in getting two closest excited states (^2P1/2,^2P3/2) of the Ir(16+) cation and their CCSD(T) correlation energies.
::

 Ir(16+) : [Xe] 4f(14) 5s(1) 5p1/2(0)  5p3/2(0)  (^2S)
 Ir(16+) : [Xe] 4f(14) 5s(0) 5p1/2(1)  5p3/2(0)  (^2P1/2)
 Ir(16+) : [Xe] 4f(14) 5s(0) 5p1/2(0)  5p3/2(1)  (^2P3/2)


For simplicity, we will work with the .X2C Hamiltonian and with the smallest Dyall.v2z decontracted basis set.

There are two ways to achieve, (i) from the converged SCF state, and, (ii) only at the correlated level from the 2P_aver SCF state.

^2P1/2,^2P3/2 SCF states
------------------------

Thanks to the symmetry with two irreps, one can place the unpaired electron either the 1st irrep (2^S1/2) or into 2nd irrep (2P1/2).

How to obtain the ^2P3/2 excited state ? For that, we use the REORDER MO keyword and read MO from the previous step.





^2P1/2,^2P3/2 CCSD(T) states
----------------------------

The other option is to start from the 2P average single determinant state in the linear symmetry, and distinguish 
between ^2P1/2,^2P3/2 states at the Coupled Cluster correlated scheme. 


