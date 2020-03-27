#!/bin/bash

# script to run tests on summit
echo "Don't forget to keep ~/.diracrc up to date"

# path for dirac.x and exacorr.x
export DIRAC_PATH=/gpfs/alpine/world-shared/chp109/DIRAC_INSTALL

# number of mpi processes (2*nnodes)
# set in pbs

if [ $# -eq 0 ]
then
  echo "give the number of the test you want to run"
else

if [ "$1" -eq 1 ] || [ "$1" -gt 8 ]
then 
  # run ccd test
  cp job_pam.pbs ccd
  export INP_file=cc_uf6.inp
  export MOL_file=UF6.xyz
  cd ccd
    bsub job_pam.pbs 
  cd ..
  echo "ccd submited"
fi

if [ "$1" -eq 2 ] || [ "$1" -gt 8 ]
then
  # run ccsd test
  cp job_pam.pbs ccsd
  export INP_file=cc_uf6.inp
  export MOL_file=UF6.xyz  
  cd ccsd
    bsub job_pam.pbs 
  cd ..
  echo "ccsd submited"
fi

if [ "$1" -eq 3 ] || [ "$1" -gt 8 ]
then
  # run scheme test
  cp job_pam.pbs scheme
  export INP_file=cc_uf6.inp
  export MOL_file=UF6.xyz
  cd scheme
    bsub job_pam.pbs
  cd ..
  echo "scheme submited"
fi

if [ "$1" -eq 4 ] || [ "$1" -gt 8 ]
then
  # run prop test
  cp job_pam.pbs prop
  export INP_file=cc_uf6.inp
  export MOL_file=UF6.xyz
  cd prop
    bsub job_pam.pbs
  cd ..
  echo "prop submited"
fi

fi
