Hartree-Fock and DFT calculations in solution with the polarizable continuum model
==================================================================================

This tutorial will cover the set up of an Hartree-Fock and a DFT calculation using the
polarizable continuum model (PCM) to include solvent effects. The implementation
of the PCM-SCF algorithm in DIRAC is documented in the recent paper [DiRemigio2014]_,
to which the reader is referred for further technical details.
The current implementation makes use of an external module for the PCM-related
tasks. The module needs some input information to run properly.

As explained in the star(PCMSOL) section of the reference manual, a PCM calculation can be set up in two ways:
* using **three** input files. The additional input file contains directives for the PCMSolver module;
* using **two** input files. The additional directives in the inp file are forwarded to the PCMSolver module.
The second method is not as flexible as the first one and we will cover both in this tutorial.

For further details, please refer to the PCM documentation in the reference manual.


An Hartree-Fock calculation
---------------------------


