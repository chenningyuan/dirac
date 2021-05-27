#!/bin/sh
#
# Script for NEVERCODE CI-service
# see  https://developer.nevercode.io/docs/custom-build-steps
#
# Add "NC_POST_CLONE_SCRIPT=./maintenance/cdash/nevercode.sh" 
# into your environment setting
#


# Show system info
echo "System info, lscpu :"
lscpu
echo "Disk, df -h :"
df -h
echo "Memory, free -m :"
free -g -h

echo "git whereabouts:"
which git
whereis git
git --version
echo

# install packages
echo "installing packages :"
#sudo apt-get update
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

# install most recent CMake
wget --no-check-certificate https://cmake.org/files/v3.5/cmake-3.5.2-Linux-x86_64.tar.gz
tar xzf cmake-3.5.2-Linux-x86_64.tar.gz
export PATH=$PWD/cmake-3.5.2-Linux-x86_64/bin:$PATH
cmake --version


# CDash project name 
  export CTEST_PROJECT_NAME="DIRACext"
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
python ./setup --fc=mpif90 --cc=mpicc --cxx=mpicxx --int64 --mpi --blas=off --lapack=off --cmake-options="-D BUILDNAME='mpich-gnu-i8_nevercode' -D DART_TESTING_TIMEOUT=9999 -D ENABLE_BUILTIN_BLAS=ON -D ENABLE_BUILTIN_LAPACK=ON -D ENABLE_UNIT_TESTS=ON"  build
cd build
ctest -D ExperimentalConfigure
ctest -j $N_PROC -D ExperimentalBuild
echo "$(ctest -j $TEST_N_PROC -D ExperimentalTest -L short)"
ctest -D ExperimentalSubmit 

# finally do the quick documentation buildup in the build_semaphore-ci_i8 directory
#make html
#make doxygen
#make slides

