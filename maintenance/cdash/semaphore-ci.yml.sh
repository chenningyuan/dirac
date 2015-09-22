#!/usr/bin/env bash
#
# Iest sript for CI build based on https://semaphoreci.com
#
# It is recommended to copy (and modify) commands from this file to "Test Settings" of your project on semaphore.
# You can also run this script when you add to your projects "Test Settings" this command: 
#
#           "chmod a+x maintenance/cdash/semaphore-ci.yml.sh && ./maintenance/cdash/semaphore-ci.yml.sh".
#
# But this is not good for debugging because everything is run as one command.
#
# **Note**
# In the following command: 
#
# echo "$(ctest -j 8 --track Miro -D ExperimentalTest -L short)"
#
# the command "echo" makes possible to continue after one of Dirac tests fails
# because otherwise server stops and do not submit to dashboard.
#


# Show system info
echo "System info, lscpu :"
lscpu
echo "Disk, df -h :"
df -h
echo "Memory, free -m :"
free -m

echo "git whereabouts:"
which git
whereis git
git --version
echo

# install packages
echo "installing packages :"
sudo apt-get update
sudo apt-get install gcc -y
sudo apt-get install g++ -y
sudo apt-get install gfortran -y
sudo apt-get install make -y
sudo apt-get install cmake -y
sudo apt-get install libblas-dev -y
sudo apt-get install liblapack-dev -y
sudo apt-get install zlib1g-dev -y
sudo apt-get install doxygen -y
sudo apt-get install python-dev -y
sudo apt-get install python-pip -y
# take care of MPICH packages
sudo apt-get remove openmpi-bin libopenmpi-dev -y
sudo apt-get install mpich2 -y
#
sudo pip install matplotlib
sudo pip install sphinx==1.2.3
sudo pip install sphinxcontrib-bibtex
sudo pip install sphinx_rtd_theme
sudo pip install hieroglyph

# set for parallel buildup
  DIRAC_MPI_COMMAND="mpirun -np 2"
# determine number of processors
  N_PROC=$(cat /proc/cpuinfo | grep processor | wc -l)
# number of available processors/jobs for test will be 2 times less because we use "mpirun -np 2"
  TEST_N_PROC=$(expr $(cat /proc/cpuinfo | grep processor | wc -l) / 2)
# print
echo "Number of processors/jobs for build step is $N_PROC"
echo "Number of processors/jobs for test step is $TEST_N_PROC"

# update submodules
git submodule update --init --recursive

# 64-bit Integer, MPICH parallel build
python ./pamadm -a --builddir=build_semaphore-ci_i8
python ./setup.py --fc=mpif90 --cc=mpicc --cxx=mpicxx --int64 --mpi --blas=off --lapack=off --cmake-options="-D BUILDNAME='mpich-gnu-i8_semaphore' -D DART_TESTING_TIMEOUT=9999 -D ENABLE_BUILTIN_BLAS=ON -D ENABLE_BUILTIN_LAPACK=ON"  build_semaphore-ci_i8
cd build_semaphore-ci_i8
ctest --track master -D ExperimentalConfigure
ctest -j $N_PROC --track master -D ExperimentalBuild
#echo "$(ctest -j 8 --track Miro -D ExperimentalTest -L short)"
echo "$(ctest -j $TEST_N_PROC --track master -D ExperimentalTest)"
ctest --track master -D ExperimentalSubmit 

# finally do the quick documentation buildup in the build_semaphore-ci_i8 directory
make html
make doxygen
make slides

