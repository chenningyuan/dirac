@echo off

rem Testing script for Cygwin which can be (directly) run from Windows Task Scheduler
rem /home/milias/test directory must exist before script is run
rem /home/milias/cdash-testing.log file contains log from previous run
rem this script can be placed everywhere you want and does not need any arguments

rem unset PUTTYs environment variable (it makes problems to Cygwin git)
set GIT_SSH=

rem delete old /home/milias/cdash-testing.log
C:\cygwin64\bin\sh --login -c "rm -f /home/milias/cdash-testing.log"

rem delete previously downloaded sources and old builds
C:\cygwin64\bin\sh --login -c "(echo -e \" deleting sources..\" && rm -rf /home/milias/test/trunk_cloned) 2>&1 | tee -a /home/milias/cdash-testing.log"

rem clone sources from gitlab
C:\cygwin64\bin\sh --login -c "(echo -e \"\n cloning sources..\" && cd /home/milias/test && git clone git@gitlab.com:dirac/dirac.git trunk_cloned) 2>&1 | tee -a /home/milias/cdash-testing.log"

rem update submodules
C:\cygwin64\bin\sh --login -c "(echo -e \"\n updating submodules..\" && cd /home/milias/test/trunk_cloned && git submodule update --init) 2>&1 | tee -a /home/milias/cdash-testing.log" 

rem configure build : 64-bit with Open MPI
rem also define ZLIB_ROOT to find correct zlib library
C:\cygwin64\bin\sh --login -c "(echo -e \"\n running setup..\" && cd /home/milias/test/trunk_cloned && python setup -D BUILDNAME="WinS12-cygwin-i8-openmpi" --int64 --mpi --cc=mpicc --cxx=mpicxx --fc=mpifort -D ZLIB_ROOT=/usr/lib build) 2>&1 | tee -a /home/milias/cdash-testing.log"
C:\cygwin64\bin\sh --login -c "(echo -e \"\n ExperimentalConfigure..\" &&cd /home/milias/test/trunk_cloned/build && ctest -D ExperimentalConfigure) 2>&1 | tee -a /home/milias/cdash-testing.log"

rem build
C:\cygwin64\bin\sh --login -c "(echo -e \"\n ExperimentalBuild..\" && cd /home/milias/test/trunk_cloned/build && ctest -j $(cat /proc/cpuinfo | grep processor | wc -l) -D ExperimentalBuild) 2>&1 | tee -a /home/milias/cdash-testing.log"

rem export DIRAC_MPI_COMMAND and test
C:\cygwin64\bin\sh --login -c "(echo -e \"\n ExperimentalTest..\" && export DIRAC_MPI_COMMAND=\"mpirun -np 2\" && cd /home/milias/test/trunk_cloned/build && ctest -j $(expr $(cat /proc/cpuinfo | grep processor | wc -l) / 4) -D ExperimentalTest) 2>&1 | tee -a /home/milias/cdash-testing.log"

rem submit to cdash
C:\cygwin64\bin\sh --login -c "(echo -e \"\n ExperimentalSubmit..\" &&cd /home/milias/test/trunk_cloned/build && ctest -D ExperimentalSubmit) 2>&1 | tee -a /home/milias/cdash-testing.log"

rem wait for user input to close this terminal window
rem pause
