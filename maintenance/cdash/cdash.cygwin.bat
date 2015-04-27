@echo off

rem Testing script for Cygwin which can be run from Windows Task Scheduler

rem delete previously downloaded sources and old builds
C:\cygwin64\bin\sh --login -c "rm -rf /home/milias/test/trunk_cloned"

rem clone sources from gitlab
C:\cygwin64\bin\sh --login -c "cd /home/milias/test && git clone git@gitlab.com:dirac/dirac.git trunk_cloned"

rem update submodules
C:\cygwin64\bin\sh --login -c "cd /home/milias/test/trunk_cloned && git submodule update --init" 

rem configure build : 64-bit with Open MPI
C:\cygwin64\bin\sh --login -c "cd /home/milias/test/trunk_cloned && python setup --int64 --mpi --cc=mpicc --cxx=mpicxx --fc=mpifort build"
C:\cygwin64\bin\sh --login -c "cd /home/milias/test/trunk_cloned/build && ctest -D ExperimentalConfigure"

rem build
C:\cygwin64\bin\sh --login -c "cd /home/milias/test/trunk_cloned/build && ctest -j $(cat /proc/cpuinfo | grep processor | wc -l) -D ExperimentalBuild"

rem export DIRAC_MPI_COMMAND and test
C:\cygwin64\bin\sh --login -c "export DIRAC_MPI_COMMAND=\"mpirun -np 2\" && cd /home/milias/test/trunk_cloned/build && ctest -j $(expr $(cat /proc/cpuinfo | grep processor | wc -l) / 4) -D ExperimentalTest"

rem submit to cdash
C:\cygwin64\bin\sh --login -c "cd /home/milias/test/trunk_cloned/build && ctest -D ExperimentalSubmit"

