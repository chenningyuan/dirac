#!/usr/bin/env bash

# Sequential test script for https://codeship.com/

# It is recommended to copy (and modify) commands from this file to "Test Settings" of your project on codeship.
# You can also run this script when you add to your projects "Test Settings" this command: 
#
#                "chmod a+x maintenance/cdash/codeship.yml.sh && ./maintenance/cdash/codeship.yml.sh".
#
# But this is not good for debugging because everything is run as one command.
#
# **Important**
# If you use to run this script command "chmod a+x maintenance/cdash/codeship.yml.sh && ./maintenance/cdash/codeship.yml.sh" it is possible 
# that even in case some command fails your build on codeship might be green - in success state!
#
# **Note**
# In the following command: 
#
#           echo "$(ctest -j $N_PROC -D ExperimentalTest -L short)"
#
# command "echo" makes possible to continue after one of Dirac tests fails because otherwise 
# codeship can stop and do not submit results to dashboard.
#
# **Note**
# https://codeship.com/documentation/faq/limit_builds/
# "We build all branches you push to your repository. In our opinion every push to your repository should be tested.
# We don't have a feature to limit which branches can be built."
#
# But it is possible to ignore or run specific command on branch: 
# see https://codeship.com/documentation/faq/ignore-command-on-branch/
#
# Written by Ivan Hrasko, Miroslav Ilias
#
#

# export number of processors as environment variable N_PROC
export N_PROC=$(cat /proc/cpuinfo | grep processor | wc -l)

# export CDash project name
export CTEST_PROJECT_NAME=DIRACext

# show number of processors to be used
echo "Number of available processors is $N_PROC"

# show system info (processors, disks, memory)
lscpu
df -h
free -m

# **Note**
# It is not possible to install packages with "sudo apt-get install".
# In case some package is missing you can write a message to codeship maintainers.
# install python packages (must be without sudo)
# note: these packages are probably already installed
pip install matplotlib
pip install sphinx==1.2.3
pip install sphinxcontrib-bibtex
pip install sphinx_rtd_theme
pip install hieroglyph

# install most recent CMake
wget --no-check-certificate https://cmake.org/files/v3.5/cmake-3.5.2-Linux-x86_64.tar.gz
tar xzf cmake-3.5.2-Linux-x86_64.tar.gz
export PATH=$PWD/cmake-3.5.2-Linux-x86_64/bin:$PATH
cmake --version


# git whereabouts
which git
whereis git
git --version

# documentation generators whereabouts
sphinx-build --version
doxygen --version

# update submodule
git submodule update --init --recursive

# Do the 32-bit Integer build
python ./setup --cmake-options="-D BUILDNAME='gnu.i4-codeship-ci' -D ENABLE_BUILTIN_LAPACK=ON -D ENABLE_UNIT_TESTS=ON" build_codeship_i4
cd build_codeship_i4
ctest -D ExperimentalConfigure
ctest -j $N_PROC -D ExperimentalBuild
echo "$(ctest -j $N_PROC  -D ExperimentalTest -L short)"
#echo "$(ctest -j $N_PROC -D ExperimentalTest)"
if [ "$CI_BRANCH" == "master" ]; then ctest -D ExperimentalSubmit; fi

# Back to the source dir
cd .. 

# Now do the 64-bit Integer build
python ./setup --int64 --blas=off --lapack=off --cmake-options="-D BUILDNAME='gnu.i8-codeship-ci' -D ENABLE_BUILTIN_BLAS=ON -D ENABLE_BUILTIN_LAPACK=ON -D ENABLE_UNIT_TESTS=ON" build_codeship_i8
cd build_codeship_i8
ctest -D ExperimentalConfigure
ctest -j $N_PROC -D ExperimentalBuild

echo "$(ctest -j $N_PROC -D ExperimentalTest -L short)"
#echo "$(ctest -j $N_PROC -D ExperimentalTest)"

if [ "$CI_BRANCH" == "master" ]; then ctest -D ExperimentalSubmit; fi

# finally do the documentation buildup in the build_codeship_i8 directory
make html
make slides
# for doxygen encompass the building command into echo so that the build does not stop in the case of error
echo "$(make doxygen)"

cd ..
