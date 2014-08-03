=========================================
The :math:`UF_6` anion - another approach
=========================================

Different and more robust solution for obtaining the SCF convergence of the :math:`UF_6` anion is provided here.

It is based on the damping of the Fock matrix from iteration to iteration.

In the first step we follow the previous approach of starting from the closed-shell UF6 system.
We also used a shift in virtual space to reduce mixing between occupied and virtual space.

We compare two calculations starting from neutral molecular spinors of UF6. 
One calculation uses the DIIS method and the other uses the damping with the damping factor of 0.90.

The damping scheme uses more of the input Fock matrix that the output Fock matrix.

Concerning the damping factor, the lower than 1, we shoudl get faster convergence. When it is getting too small,
iterations diverge. The closer we are to the converged states, the smalle damping factor can be.

With damping convergence speedup the calculation converges to 1,4 10^-8 after 50 iterations.
On the other hand the DIIS iterations initially converge to about 5*10^-5 but tthen began to diverge
and after 150 iterations the error gradient of the Fock matrix is about 2.10^-5 .

%On the other hand if we began from the converged solutions from the damping, the DIIS iterations converge 
%rapidly after 7 iterations to 8.10^-9

TODO: automatic adjustement o the damping factor, suggestions for keyword
.CNVDMP
1.D-4 1.D-5 1.0D-6
.DAMPDEL ! decrease level of damp. factor
0.2



