#!/bin/bash

# example recipe to build and run tests for CDash
# you can either run this script interactively with "make Experimental" (to calibrate it)
# or execute it as a cron job with "make Nightly"
# make sure cmake is installed before you run it

# miro: modified to local update&build

# the directory where your DIRAC sits
DIRAC=/home/ilias/QCH_Work/qch_softs/dirac_space/svn_operated/g95_gcc_noopt_ownlib_ILP64/adm/trunk

# go to the local DIRAC directory
cd $DIRAC

# svn update, pamadm 
svn up 
./pamadm -a

# set environment (you may or may not need it for your compiler)
# export PATH=$PATH':/home/bast/openmpi/bin'
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH':/home/bast/openmpi/lib'
# source /opt/intel/bin/compilervars.sh ia32
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH':/opt/intel/mkl/lib/ia32'
source ~/.bashrc

# configure step
# ./setup --fc=gfortran --cc=gcc --debug --buildname="linux.gfortran.gcc.debug"
# ./setup --fc=gfortran --cc=gcc --buildname="linux.gfortran.gcc"
#  ./setup --fc=ifort --cc=icc --math="-L/opt/intel/mkl/lib/ia32 -lmkl_lapack95 -lmkl_intel -lmkl_intel_thread -lmkl_core -liomp5 -lpthread" --buildname="linux.intel12.mkl"
# ./setup --fc=mpif90 --cc=mpicc --mpi --math="-L/usr/lib -lblas -L/usr/lib -llapack -lcblas -lf77blas -latlas -L/usr/lib -llapack" --buildname="linux.gfortran.atlas.mpi"
./setup --fc=/usr/bin/g95 --cc=/usr/bin/gcc  --int64  --nomath --buildname="linux.g95.gcc.ilp64.ownlibs"

cd $DIRAC/cbuild
cmake ..

# pam expects dirac.x to be in the same directory
# you cannot copy it because it does not exist yet
# make below will create it but same step launches pam (via runtest)
ln -s $DIRAC/cbuild/dirac.x $DIRAC/dirac.x

# build and run tests
# use either Experimental or Nightly
  make Experimental -j 2
# make Nightly -j 2

exit 0
