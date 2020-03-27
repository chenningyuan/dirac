#!/bin/bash

# script to run tests on summit
echo "Don't forget to keep ~/.diracrc up to date"

# path for dirac.x and exacorr.x
export DIRAC_PATH=/gpfs/alpine/world-shared/chp109/DIRAC_INSTALL

# number of mpi processes (2*nnodes)
export N_NODES=30

# run ccd test
cp job_pam.pbs ccd
export INP_file=UF6.xyz
export MOL_file=cc_uf6.inp
cd ccd
  bsub job_pam.pbs 
cd ..
echo "ccd submited"

