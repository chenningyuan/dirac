====
MnO6 
====

The MnO6 compound is good example of the oxide crystal structure.

In the following we show how one can get converged molecular spinors.
We take into account two parameters for controling the convergence -
the OPENFAC and OLEVEL.

The noAcAc abbreviations means OPENFAC=0, what is no active-active interation among open-shell
spinors.

The wAcAc means OPENFAC=1, that is full  active-active interation (or full open-shell contribution 
to the Fock matrix).

1. Starting from NoAcAc we get convergence. We save the DFCOEF file with molecular spinors.

2. Using AcAc with molecular spinors (DFCOEF file) from previous SCF run does not converge.

3. Using AcAc with OLEVEL=0.2 and with DFCOEF from point 1 does converge.

