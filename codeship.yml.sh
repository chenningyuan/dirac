#!/usr/bin/env bash

# Sequential test script for https://codeship.com/

# It is recommended to copy (and modify) commands from this file to "Test Settings" of your project on codeship.
# You can also run this script when you add to your projects "Test Settings" this command: "chmod a+x codeship.yml.sh && ./codeship.yml.sh".
# But this is not good for debugging because everything is run as one command.

# **Important**
# If you use to run this script command "chmod a+x codeship.yml.sh && ./codeship.yml.sh" it is possible that even in case some command fails your build
# on codeship might be green - in success state!

# **Note**
# In the following command: 
# echo "$(ctest -j $N_PROC --track Miro -D ExperimentalTest -L short)"
# command "echo" makes possible to continue after one of Dirac tests fails because otherwise codeship can stop and do not submit to dashboard.

# **Note**
# https://codeship.com/documentation/faq/limit_builds/
# "We build all branches you push to your repository. In our opinion every push to your repository should be tested.
# We don't have a feature to limit which branches can be built."
#
# But it is possible to ignore or run specific command on branch: 
# see https://codeship.com/documentation/faq/ignore-command-on-branch/

# export number of processors as environment variable N_PROC
export N_PROC=$(cat /proc/cpuinfo | grep processor | wc -l)

# show number of processors to be used
echo "Number of available processors is $N_PROC"

# Show system info
lscpu
df -h
free -m

# git whereabouts
which git
whereis git
git --version

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

# 32-bit Integer build
python ./pamadm -a --builddir=build_codeship_i4
python ./setup -D BUILDNAME="gnu.i4-codeship-ci" build_codeship_i4
cd build_codeship_i4
ctest --track Miro -D ExperimentalConfigure
ctest -j $N_PROC --track Miro -D ExperimentalBuild
echo "$(ctest -j $N_PROC --track Miro -D ExperimentalTest -L short)"
if [ "$CI_BRANCH" == "master" ]; then ctest --track Miro -D ExperimentalSubmit; fi

# Back to source dir
cd .. 

# 64-bit Integer build
python ./pamadm -a --builddir=build_codeship_i8
python ./setup -D BUILDNAME="gnu.i8-codeship-ci" --int64 build_codeship_i8
cd build_codeship_i8
ctest --track Miro -D ExperimentalConfigure
ctest -j $N_PROC --track Miro -D ExperimentalBuild
echo "$(ctest -j $N_PROC --track Miro -D ExperimentalTest -L short)"
if [ "$CI_BRANCH" == "master" ]; then ctest --track Miro -D ExperimentalSubmit; fi

# finally do the documentation buildup in the build_codeship_i8 directory
make html
make doxygen
make slides
