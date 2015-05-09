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

# show system info (processors, disks, memory)
lscpu
df -h
free -m
# show version info for Open MPI
ompi_info

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

# git whereabouts
which git
whereis git
git --version

# documentation generators whereabouts
sphinx-build --version
doxygen --version

# configure Dirac with Open MPI in 64-bit mode
python ./pamadm -a --builddir=build
python ./setup -D BUILDNAME="open_mpi.gnu.i8-codeship-ci" --int64 --mpi -D DART_TESTING_TIMEOUT=99999 build
cd build
ctest -D ExperimentalConfigure
# build
ctest -j $(cat /proc/cpuinfo | grep processor | wc -l) -D ExperimentalBuild
# show linked libraries
ldd dirac.x
# run ctest (note: (ctest ..) || echo "CTest end." makes possible to continue after some tests inside ctest command fail)
(ctest -j $(expr $(cat /proc/cpuinfo | grep processor | wc -l) / 4) -D ExperimentalTest) || echo "CTest end."
# submit
ctest -D ExperimentalSubmit

# finally build the documentation
make html
make slides
# for doxygen encompass the building command into echo so that the build does not stop in the case of error
# **Note**
# Using echo command like in this case is not enough for ctest because it can happen that command prints no output
# for long time and build is stopped.
echo "$(make doxygen)"
