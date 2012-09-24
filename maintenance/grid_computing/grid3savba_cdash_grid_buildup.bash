#!/bin/bash
#############################################################################
#
#      parametrized script for DIRAC buildup in the grid environment
#
#############################################################################

# check parameter - VO name

if [[ $1 != "voce" && $1 != "compchem" && $1 != "isragrid" && $1 != "osg" ]]; then
 echo -e "\n wrong parameter - VO name : $1 "
 exit -12
else
 VO=$1
 echo -e "OK, you specified properly VO=$VO"
 exit 0
fi

echo -e "\n *** Working host is: "; hostname -f

source $HOME/intel/bin/compilervars.sh intel64
echo "Intel Fortran/C/C++ noncommercial compilers with MKL library activated."
# My own OpenMPI Intel static
export PATH=/home/ilias/bin/openmpi_i32lp64_intel_static/bin:$PATH
export LD_LIBRARY_PATH=/home/ilias/bin/openmpi_i32lp64_intel_static/lib:$LD_LIBRARY_PATH
# my cmake stuff
export PATH=/home/ilias/bin/cmake_static:$PATH
echo "My own CMake is in PATH."
# my most recent python 2.7.1 #
#export PATH=/people/disk2/ilias/bin/python/Python-2.7.1:$PATH
echo "My PATH=$PATH"
#echo "which make:"; which make

timestamp1=`date +\%F_\%k-\%M-\%S`; timestamp=${timestamp1// /};
echo "Running DIRAC cdash buildup at "$timestamp

DIRAC=/home/ilias/qch_work/software/dirac_git/trunk_current
#
BUILD=build_Intel_OpenMPI_MKL_ILP64_STATIC
BUILD_SERIAL=build_Intel_serial_MKL_ILP64_STATIC
BUILD_GNU_SERIAL=build_GNU_serial_ILP64_STATIC

cd $DIRAC

export MATH_ROOT=$HOME/intel/mkl; echo "Activated MATH_ROOT=$MATH_ROOT"

# download fresh clone...
git pull

./setup --fc=/home/ilias/bin/openmpi_i32lp64_intel_static/bin/mpif90 --cc=/home/ilias/bin/openmpi_i32lp64_intel_static/bin/mpicc --cxx=/home/ilias/bin/openmpi_i32lp64_intel_static/bin/mpiCC --static --mpi -D BUILDNAME="grid3_savba_ompi_intel_mkl_ilp64_STATIC" -D RUNTEST_NR_PROC=2 --int64 -D RUNTEST_TIMEOUT="00:18:00" -D DART_TESTING_TIMEOUT=99999 -D VERBOSE_OUTPUT=ON -D ENABLE_MAGNETIZ_FIX=OFF -D RUNTEST_BINARY=$DIRAC/$BUILD/dirac.x --force $BUILD

cd $DIRAC/$BUILD
 #cmake ..
 #make -j2 dirac.x
 make depend
 #make ExperimentalUpdate
 make ExperimentalConfigure
 make ExperimentalBuild
 make ExperimentalSubmit

 ls -lt $DIRAC/$BUILD/dirac.x
 ldd $DIRAC/$BUILD/dirac.x

timestamp1=`date +\%F_\%k-\%M-\%S`; timestamp=${timestamp1// /};
echo -e "\n build 1 - Intel-OpenMPI,MKL,static - finished at $timestamp \n"

cd $DIRAC
./setup --fc=ifort --cc=icc --cxx=icpc --static -D BUILDNAME="grid3_savba_serial_intel_mkl_ilp64_STATIC" --int64 -D RUNTEST_TIMEOUT="00:18:00" -D DART_TESTING_TIMEOUT=99999 -D VERBOSE_OUTPUT=ON -D ENABLE_MAGNETIZ_FIX=ON -D RUNTEST_BINARY=$DIRAC/$BUILD_SERIAL/dirac.x --force $BUILD_SERIAL
cd $DIRAC/$BUILD_SERIAL
 #cmake ..
 #make -j2 dirac.x
 make depend
 #make ExperimentalUpdate
 make ExperimentalConfigure
 make ExperimentalBuild
 make ExperimentalSubmit

 ls -lt $DIRAC/$BUILD_SERIAL/dirac.x
 ldd $DIRAC/$BUILD_SERIAL/dirac.x

timestamp1=`date +\%F_\%k-\%M-\%S`; timestamp=${timestamp1// /};
echo -e "\n build 2 - Intel,MKL,serial,static - finished at $timestamp \n"

# cancel the variable #
unset MATH_ROOT
cd $DIRAC
./setup --fc=gfortran44 --cc=gcc44 --cxx=g++44 --static -D BUILDNAME="grid3_savba_serial_GNU_ilp64_STATIC" --int64 -D RUNTEST_TIMEOUT="00:18:00" -D DART_TESTING_TIMEOUT=99999 -D VERBOSE_OUTPUT=ON -D ENABLE_MAGNETIZ_FIX=OFF -D RUNTEST_BINARY=$DIRAC/$BUILD_GNU_SERIAL/dirac.x --force $BUILD_GNU_SERIAL
cd $DIRAC/$BUILD_GNU_SERIAL
 #cmake ..
 #make -j2 dirac.x
 make depend
 #make ExperimentalUpdate
 make ExperimentalConfigure
 make ExperimentalBuild
 make ExperimentalSubmit

 ls -lt $DIRAC/$BUILD_GNU_SERIAL/dirac.x
 ldd $DIRAC/$BUILD_GNU_SERIAL/dirac.x
timestamp1=`date +\%F_\%k-\%M-\%S`; timestamp=${timestamp1// /};
echo -e "\n build 3 - GNU,serial,static -  finished at $timestamp \n"

#################################################
### continue with the launching for the grid  ###
#################################################

# create the reduced tar-file containing all important stuff ...
cd $DIRAC
# copy there the cmake static
cp -R /home/ilias/bin/cmake_static/*                               $DIRAC/$BUILD/.
cp -R /home/ilias/bin/cmake_static/*                               $DIRAC/$BUILD_SERIAL/.

# copy there the mpirun static with related files
cp    /home/ilias/bin/openmpi_i32lp64_intel_static/bin/mpirun      $DIRAC/$BUILD/.
/bin/rm -rf $BUILD/share
mkdir $BUILD/share
cp -R -P /home/ilias/bin/openmpi_i32lp64_intel_static/share/openmpi   $DIRAC/$BUILD/share/.

# problems with "--exclude test/*/*.out"  - so remove them as such
/bin/rm test/*/*.out*
/bin/rm -r test/*/compare

echo -e "\n now packing stuff for grid computations..."

tar czf dirac_grid_suite.tgz   CMakeLists.txt CTestConfig.cmake pamadm  pam runtest setup  test  basis  basis_dalton basis_ecp  cmake $BUILD $BUILD_SERIAL doc src utils COPYING maintenance/grid_computing/diracrc_grid

echo -e "\n packing done, we have tar-ball file ready to be put on grid's SE : "
ls -lt dirac_grid_suite.tgz
#echo -e "\n This tarball contains :"
#tar -tvf dirac_grid_suite.tgz

######################################################## 
#
# initialize your certificate from grid computing
#
# voms-proxy-init --voms voce      -vomslife 24:00 
# voms-proxy-init --voms compchem  -vomslife 24:00 
# 
######################################################## 

echo -e "\n ---  My valid certificate (voms-proxy-info --all)  ---"
voms-proxy-info --all ||  { echo -e "\n Exiting...no valid certificate found !\n" ; exit -2; }

# delete the old package from the SE
echo -e "\n deleting the previous file, lcg-del -a -v lfn:/grid/$VO/ilias/dirac_grid_suite.tgz ...."
lcg-del -a -v  lfn:/grid/$VO/ilias/dirac_grid_suite.tgz ||  { echo "Exit, deleting failed!" ; exit -3; }
echo -e "\n After deleting, checking emptyness of VO's storage space; lcg-ls -l lfn://grid/$VO/ilias"
lcg-ls -l lfn://grid/$VO/ilias || { echo "Exit, lcg-ls command failed !" ; exit -4; }

# upload the fresh tar-ball onto the SE
echo -e "\n uploading lcg-cr -v file:$PWD/dirac_grid_suite.tgz -l lfn://grid/voce/ilias/dirac_grid_suite.tgz...."
lcg-cr -v file:$PWD/dirac_grid_suite.tgz  -l lfn://grid/$VO/ilias/dirac_grid_suite.tgz || { echo "Exit, uploading failed !" ; exit -4; }
#lcg-cr file:$PWD/dirac_grid_suite.tgz  -l lfn://grid/compchem/ilias/dirac_grid_suite.tgz

###   Check your copied file on the grid SE:  ###
echo -e "\n After uploading, list of files on grid SE, lcg-ls -l lfn://grid/$VO/ilias/dirac_grid_suite.tgz"
lcg-ls -l lfn://grid/$VO/ilias/dirac_grid_suite.tgz
#  lcg-ls lfn://grid/compchem/ilias
echo -e "\n File attributes lcg-lr  lfn://grid/$VO/ilias/dirac_grid_suite.tgz:"
lcg-lr    lfn://grid/$VO/ilias/dirac_grid_suite.tgz
#  do also uberftp se.ui.savba.sk "ls /dpm/ui.savba.sk/home...."

#############################################################################
#            submit external launching script with glite
#############################################################################
timestamp1=`date +\%F_\%k-\%M-\%S`; timestamp=${timestamp1// /};
echo -e "\n going to launch for grid, with timestamp=$timestamp..."
cd maintenance/grid_computing
echo -e "I am in the directory $PWD"

glite-wms-job-submit -o JOB.$timestamp.$VO.1 -a submit_$VO.jdl
glite-wms-job-submit -o JOB.$timestamp.$VO.2 -a submit_$VO.jdl
glite-wms-job-submit -o JOB.$timestamp.$VO.3 -a submit_$VO.jdl

echo -e "\n Done, script finished."
exit 0
