
**************************
Plotting orbital densities
**************************

In this tutorial we will try to plot many orbital densities at once in a nice
and lazy way. We will use the script wrapper.py to achieve this.

To keep things simple and fast we will calculate the Ne atom in a cc-pVDZ basis.
But the principle is the same for (large) molecules.

Getting the wave function
=========================

Let us use the following molecule input (ne.mol)::

  INTGRL


  C   1
         10.    1
  Ne     0.0  0.0  0.0
  LARGE BASIS cc-pVDZ
  FINISH

Together with a simple job input (scf.inp)::

  **DIRAC
  .WAVE FUNCTION
  **HAMILTONIAN
  .DFT
   LDA
  **WAVE FUNCTION
  .SCF
  *END OF INPUT

We can get the wave function with the following run script::

  #!/bin/bash

  # adapt this
  DIRAC_ROOT=$HOME/dirac
  export BASDIR=$DIRAC_ROOT/basis_dalton
  SCRATCH=$HOME/trash/scratch

  # define command
  WRAPPER=$DIRAC_ROOT/utils/wrapper.py
  EXECUTABLE=$DIRAC_ROOT/build/dirac.x
  dirac="python $WRAPPER --executable=$EXECUTABLE --tmpfull=$SCRATCH"

  # get DFCOEF ("wave function")
  $dirac ne.mol scf.inp -gDFCOEF

  # clean up
  rm -rf $SCRATCH
  rm -f *tgz

  exit 0

After running the script you should see the file DFCOEF in your submit directory.

Plotting the orbitals
=====================

Ne (cc-pVDZ) has 8 gerade and 6 ungerade orbitals. Let us get them plotted to cube files
using the following job input (plot.inp)::

  **DIRAC
  #.WAVE FUNCTION
  **HAMILTONIAN
  .DFT
   LDA
  **WAVE FUNCTION
  .SCF
  **VISUAL
  .3D
   80 80 80         # will generate 80x80x80 cubes
  .DENSITY
   DFCOEF
  .ORBSTRING        # expects two lines: gerade and ungerade
   orb_g..orb_g
   orb_u..orb_u
  .3D_INT           # to double-check that density integrates to 2 electrons
  *END OF INPUT

And the following run script::

  #!/bin/bash

  # adapt this
  DIRAC_ROOT=$HOME/dirac
  export BASDIR=$DIRAC_ROOT/basis_dalton
  SCRATCH=$HOME/trash/scratch

  # define command
  WRAPPER=$DIRAC_ROOT/utils/wrapper.py
  EXECUTABLE=$DIRAC_ROOT/build/dirac.x
  dirac="python $WRAPPER --executable=$EXECUTABLE --tmpfull=$SCRATCH"

  # get all gerade orbitals
  orb_u=0
  for orb_g in `seq 1 8`; do
      $dirac ne.mol plot.inp -iDFCOEF -rorb_g=$orb_g -rorb_u=$orb_u \
                                      -gplot.3d.cube=${orb_g}_${orb_u}.cube
  done

  # get all ungerade orbitals
  orb_g=0
  for orb_u in `seq 1 6`; do
      $dirac ne.mol plot.inp -iDFCOEF -rorb_g=$orb_g -rorb_u=$orb_u \
                                      -gplot.3d.cube=${orb_g}_${orb_u}.cube
  done

  # clean up
  rm -rf $SCRATCH
  rm -f *tgz

  exit 0

After running this script you should see many cube files in your submit directory::

  -rw------- 1 radovan radovan 824K Aug 16 12:52 0_1.cube
  -rw------- 1 radovan radovan 824K Aug 16 12:52 0_2.cube
  -rw------- 1 radovan radovan 824K Aug 16 12:52 0_3.cube
  -rw------- 1 radovan radovan 824K Aug 16 12:52 0_4.cube
  -rw------- 1 radovan radovan 824K Aug 16 12:52 0_5.cube
  -rw------- 1 radovan radovan 824K Aug 16 12:52 0_6.cube
  -rw------- 1 radovan radovan 824K Aug 16 12:51 1_0.cube
  -rw------- 1 radovan radovan 824K Aug 16 12:51 2_0.cube
  -rw------- 1 radovan radovan 824K Aug 16 12:51 3_0.cube
  -rw------- 1 radovan radovan 824K Aug 16 12:51 4_0.cube
  -rw------- 1 radovan radovan 824K Aug 16 12:51 5_0.cube
  -rw------- 1 radovan radovan 824K Aug 16 12:51 6_0.cube
  -rw------- 1 radovan radovan 824K Aug 16 12:51 7_0.cube
  -rw------- 1 radovan radovan 824K Aug 16 12:52 8_0.cube
