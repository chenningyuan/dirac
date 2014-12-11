Hartree-Fock and DFT calculations in solution with the polarizable continuum model
==================================================================================

This tutorial will cover the set up of an Hartree-Fock and a DFT calculation using the
polarizable continuum model (PCM) to include solvent effects. The implementation
of the PCM-SCF algorithm in DIRAC is documented in the recent paper [DiRemigio2014]_,
to which the reader is referred for further technical details.
The current implementation makes use of an external module for the PCM-related
tasks. The module needs some input information to run properly.

As explained in the `*PCMSOL` section of the reference manual, a PCM calculation can be set up in two ways:
* using **three** input files. The additional input file contains directives for the PCMSolver module;
* using **two** input files. The additional directives in the inp file are forwarded to the PCMSolver module.
The second method is not as flexible as the first one and we will cover both in this tutorial.

For further details, please refer to the PCM documentation in the reference manual.

We will use the extremely simple examples contained in the test/pcm_energy directory.
They are more than enough to show a simple workflow, possible gotchas and what to look for
in a PCM output.

An Hartree-Fock calculation
---------------------------

In this example, the first input method will be used. We have to prepare three input files.
The .mol file is prepared in the exact same manner as for an in vacuo calculation.
In this case we have the methane molecule with a STO-3G basis set::

   INTGRL
   methane (a standard MOLFDIR test at that geometry)
   nonrelativistic contracted STO-3G basis
   C   2    0          
           6.    1
   C   0.0000    0.0000    0.0000
   LARGE BASIS STO-3G 
           1.    4
   H     .629889144   .629889144  .629889144
   H    -.629889144  -.629889144  .629889144
   H     .629889144  -.629889144 -.629889144
   H    -.629889144   .629889144 -.629889144
   LARGE BASIS STO-3G
   FINISH

The .inp file will feature some additional keywords::

   **DIRAC
   .WAVE FUNCTION
   **HAMILTONIAN
   .LEVY-LEBLOND
   .PCM
   *PCM
   .SEPARATE
   .PRINT
    6
   **WAVE FUNCTION
   .SCF
   *SCF
   .PRINT
   1
   *END OF

The `.PCM` keyword enables the addition of solvent terms to the Fock matrix. Details can be
controlled through directives under the `*PCM` section (to which the reader is referred for 
further details)
