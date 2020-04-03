#!/bin/bash

# script to run tests on summit
echo "Don't forget to keep ~/.diracrc up to date"

# path for dirac.x and exacorr.x
export DIRAC_PATH=/ccs/home/jpoto/caar-zero/build_gcc811_essl_cuda_release
# export DIRAC_PATH=/gpfs/alpine/world-shared/chp109/DIRAC_INSTALL

if [ $# -eq 0 ]
then
  echo "give the number of the test you want to run"
else

if [ "$1" -eq 1 ] || [ "$1" -gt 8 ]
then
  fold=scheme1
  cp job_pam.pbs ${fold}
  export INP_file=cc_uf6.inp
  export MOL_file=UF6.xyz
  cd ${fold}
    bsub job_pam.pbs
  cd ..
  echo "${fold} submited"
fi

if [ "$1" -eq 2 ] || [ "$1" -gt 8 ]
then
  fold=scheme2
  cp job_pam.pbs ${fold}
  export INP_file=cc_uf6.inp
  export MOL_file=UF6.xyz
  cd ${fold}
    bsub job_pam.pbs
  cd ..
  echo "${fold} submited"
fi

if [ "$1" -eq 3 ] || [ "$1" -gt 8 ]
then
  fold=scheme3
  cp job_pam.pbs ${fold}
  export INP_file=cc_uf6.inp
  export MOL_file=UF6.xyz
  cd ${fold}
    bsub job_pam.pbs
  cd ..
  echo "${fold} submited"
fi

if [ "$1" -eq 4 ] || [ "$1" -gt 8 ]
then
  fold=scheme4
  cp job_pam.pbs ${fold}
  export INP_file=cc_uf6.inp
  export MOL_file=UF6.xyz
  cd ${fold}
    bsub job_pam.pbs
  cd ..
  echo "${fold} submited"
fi

if [ "$1" -eq 5 ] || [ "$1" -gt 8 ]
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

fi
