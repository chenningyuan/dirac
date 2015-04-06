#!/usr/bin/env bash

# test sript for https://semaphoreci.com

# It is recommended to copy (and modify) commands from this file to "Test Settings" of your project on semaphore.
# You can also run this script when you add to your projects "Test Settings" this command: "chmod a+x semaphore-ci.yml.sh && ./semaphore-ci.yml.sh".
# But this is not good for debugging because everything is run as one command.

# **Note**
# In the following command:
# echo "$(ctest -j 8 --track Miro -D ExperimentalTest -L short)"
# command "echo" makes possible to continue after one of Dirac tests fails because otherwise server stops and do not submit to dashboard.


# Show system info
lscpu
df -h
free -m

which git
whereis git
git --version

# install packages
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
sudo apt-get install python-matplotlib -y
sudo pip install sphinx==1.2.3
sudo pip install sphinxcontrib-bibtex
sudo pip install sphinx_rtd_theme
sudo pip install hieroglyph

# 32-bit Integer build
python ./pamadm -a --builddir=build_semaphore-ci_i4
python ./setup -D BUILDNAME="gnu.i4-semaphore-ci"  build_semaphore-ci_i4
cd build_semaphore-ci_i4
ctest --track Miro -D ExperimentalUpdate
ctest --track Miro -D ExperimentalConfigure
ctest -j 8 --track Miro -D ExperimentalBuild
echo "$(ctest -j 8 --track Miro -D ExperimentalTest -L short)"
if [ "$CI_BRANCH" == "master" ]; then ctest --track Miro -D ExperimentalSubmit; fi

# Back to source dir
cd ..

# 64-bit Integer build
python ./pamadm -a --builddir=build_semaphore-ci_i8
python ./setup -D BUILDNAME="gnu.i8-semaphore-ci" --int64 build_semaphore-ci_i8
cd build_semaphore-ci_i8
ctest --track Miro -D ExperimentalUpdate
ctest --track Miro -D ExperimentalConfigure
ctest -j 8 --track Miro -D ExperimentalBuild
echo "$(ctest -j 8 --track Miro -D ExperimentalTest -L short)"
if [ "$CI_BRANCH" == "master" ]; then ctest --track Miro -D ExperimentalSubmit; fi

# finally do the documentation buildup in the build_semaphore-ci_i8
make html
make doxygen
make slides

