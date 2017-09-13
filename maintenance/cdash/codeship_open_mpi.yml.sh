#!/usr/bin/env bash

# 64-bit parallel (Open MPI) test script for https://codeship.com/

# It is recommended to copy (and modify) commands from this file to "Test Settings" of your project on codeship.
# You can also run this script when you add to your projects "Test Settings" this command: 
#
#                "chmod a+x maintenance/cdash/codeship_open_mpi.yml.sh && ./maintenance/cdash/codeship_open_mpi.yml.sh".
#
# But this is not good for debugging because everything is run as one command.
#
# **Important**
# If you use to run this script command "chmod a+x maintenance/cdash/codeship_open_mpi.yml.sh && ./maintenance/cdash/codeship_open_mpi.yml.sh" 
# it is possible that even in case some command fails your build on codeship might be green - in success state!
#
# **Note**
# https://codeship.com/documentation/faq/limit_builds/
# "We build all branches you push to your repository. In our opinion every push to your repository should be tested.
# We don't have a feature to limit which branches can be built."
#
# But it is possible to ignore or run specific command on branch: 
# see https://codeship.com/documentation/faq/ignore-command-on-branch/

# export DIRAC_MPI_COMMAND
export DIRAC_MPI_COMMAND="mpirun -np 2"

#export CTEST_PROJECT_NAME="DIRACext"
export CTEST_PROJECT_NAME="DIRAC"

# show system info (processors, disks, memory)
lscpu
df -h
free -m
# show version info for Open MPI
ompi_info

# **Note**
# It is not possible to install packages with "sudo apt-get install".
# In case some package is missing you can write a message to codeship maintainers.

# install most recent CMake
wget --no-check-certificate https://cmake.org/files/v3.5/cmake-3.5.2-Linux-x86_64.tar.gz
tar xzf cmake-3.5.2-Linux-x86_64.tar.gz
export PATH=$PWD/cmake-3.5.2-Linux-x86_64/bin:$PATH
cmake --version

# git whereabouts
which git
whereis git
git --version

#update submodules
git submodule update --init --recursive

# configure Dirac with Open MPI in 64-bit mode
python ./setup --fc=mpif90 --cc=mpicc --cxx=mpicxx --blas=off --lapack=off --int64 --mpi --cmake-options="-D BUILDNAME='openmpi.i8.codeshipCI' -D DART_TESTING_TIMEOUT=99999  -D ENABLE_BUILTIN_BLAS=ON -D ENABLE_BUILTIN_LAPACK=ON -D ENABLE_UNIT_TESTS=ON" build_openmpi_i8
cd build_openmpi_i8
ctest -D ExperimentalConfigure
# build
ctest -j $(cat /proc/cpuinfo | grep processor | wc -l) -D ExperimentalBuild
# show linked libraries
ldd dirac.x
# run ctest (note: (ctest ..) || echo "CTest end." makes possible to continue after some tests inside ctest command fail)
(ctest -j $(expr $(cat /proc/cpuinfo | grep processor | wc -l) / 2) -D ExperimentalTest) || echo "CTest end."
# submit to CDash board
ctest -D ExperimentalSubmit

