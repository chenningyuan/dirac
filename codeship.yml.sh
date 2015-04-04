#!/usr/bin/env bash

# Sequential test sript for https://codeship.com/

# It is recommended to copy (and modify) commands from this file to "Test Settings" of your project on codeship.
# You can also run this script when you add to your projects "Test Settings" this command: "chmod a+x codeship.yml.sh && ./codeship.yml.sh".
# But this is not good for debugging because everything is run as one command.

# **Note**
# In the following command: 
# echo "$(ctest -j 16 --track Miro -D ExperimentalTest -L short)"
# command "echo" makes possible to continue after one of Dirac tests fails because otherwise codeship stops and do not submit to dashboard.

# **Note**
# https://codeship.com/documentation/faq/limit_builds/
# "We build all branches you push to your repository. In our opinion every push to your repository should be tested.
# We don’t have a feature to limit which branches can be built."
#
# But it is possible to ignore or run specific command on branch: 
# see https://codeship.com/documentation/faq/ignore-command-on-branch/


# Show system info
lscpu
df -h
free -m

# 32-bit build
python ./pamadm -a --builddir=build_codeship_i4
python ./setup -D BUILDNAME="gnu.i4-codeship-ci" build_codeship_i4
cd build_codeship_i4
ctest --track Miro -D ExperimentalUpdate
ctest --track Miro -D ExperimentalConfigure
ctest -j 32 --track Miro -D ExperimentalBuild
echo "$(ctest -j 16 --track Miro -D ExperimentalTest -L short)"
if [ "$CI_BRANCH" == "master" ]; then ctest --track Miro -D ExperimentalSubmit; fi

# Back to source dir
cd .. 

# 64-bit build
python ./pamadm -a --builddir=build_codeship_i8
python ./setup -D BUILDNAME="gnu.i8-codeship-ci" --int64 build_codeship_i8
cd build_codeship_i8
ctest --track Miro -D ExperimentalUpdate
ctest --track Miro -D ExperimentalConfigure
ctest -j 32 --track Miro -D ExperimentalBuild
echo "$(ctest -j 16 --track Miro -D ExperimentalTest -L short)"
if [ "$CI_BRANCH" == "master" ]; then ctest --track Miro -D ExperimentalSubmit; fi
