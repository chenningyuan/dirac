#!/usr/bin/env bash

# Test sript for https://codeship.com/

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
ctest --track Miro -D ExperimentalSubmit

# Back to source dir
cd .. 

#64-bit build
python ./pamadm -a --builddir=build_codeship_i8
python ./setup -D BUILDNAME="gnu.i8-codeship-ci" --int64 build_codeship_i8
cd build_codeship_i8
ctest --track Miro -D ExperimentalUpdate
ctest --track Miro -D ExperimentalConfigure
ctest -j 32 --track Miro -D ExperimentalBuild
echo "$(ctest -j 16 --track Miro -D ExperimentalTest -L short)"
ctest --track Miro -D ExperimentalSubmit
