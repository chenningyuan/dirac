#!/bin/bash
#############################################################################
#
#      parametrized script for DIRAC buildup in the grid environment
#
#############################################################################

##################################################
# fill the PATH to be the same as in online mode
##################################################
# rather do the sourcing ...
source ~/.bashrc
#source ~/.bash_profile
#source /etc/bashrc

source $HOME/intel/bin/compilervars.sh intel64
echo "Intel Fortran/C/C++ noncommercial compilers with MKL library activated."
# My own OpenMPI Intel static
export PATH=/home/ilias/bin/openmpi_intel_stat/bin:$PATH
export LD_LIBRARY_PATH=/home/ilias/bin/openmpi_intel_stat/lib:$LD_LIBRARY_PATH
# my cmake stuff
export PATH=/home/ilias/bin/newer_cmake:$PATH
echo "My own (dynamic) CMake is in PATH."

echo "My PATH=$PATH"
#echo "which make:"; which make
echo -e "\n set:"; set

timestamp1=`date +\%F_\%k-\%M-\%S`; timestamp=${timestamp1// /};
echo "Running DIRAC cdash buildup at "$timestamp

DIRAC=/home/ilias/Work/qch/software/trunk
#
BUILD=build_intel_openmpi_static
BUILD_SERIAL=build_serial_intel_static


cd $DIRAC

export MATH_ROOT=$HOME/intel/mkl; echo "Activated MATH_ROOT=$MATH_ROOT"

# update the local trunk...
echo -e "\n updating the local $DIRAC"
git pull

if [[ -d "$BUILD" ]]; then
   echo "removing previous build directory $BUILD"
  /bin/rm -rf $BUILD
fi
./setup --fc=mpif90 --cc=mpicc --cxx=mpiCC --static -D BUILDNAME="grid_savba_ompi_intel_mkl_ilp64_STATIC" -D RUNTEST_NR_PROC=2 --int64 -D RUNTEST_TIMEOUT="00:18:00" -D DART_TESTING_TIMEOUT=99999  $BUILD
 cd $BUILD
 make ExperimentalUpdate
 make ExperimentalConfigure
 make ExperimentalBuild
 make ExperimentalSubmit

 ls -lt $DIRAC/$BUILD/dirac.x
 ldd $DIRAC/$BUILD/dirac.x
 python binary_info.py > VERSION_cmake  2>&1  # get version

timestamp1=`date +\%F_\%k-\%M-\%S`; timestamp=${timestamp1// /};
echo -e "\n build 1 - Intel-OpenMPI,MKL,static - finished at $timestamp \n"

cd $DIRAC
if [[ -d "$BUILD_SERIAL" ]]; then
   echo "removing previous build directory $BUILD_SERIAL"
  /bin/rm -rf $BUILD_SERIAL
fi

./setup --fc=ifort --cc=icc --cxx=icpc --static -D BUILDNAME="grid_savba_serial_intel_mkl_ilp64_STATIC" --int64 -D RUNTEST_TIMEOUT="00:18:00" -D DART_TESTING_TIMEOUT=99999  $BUILD_SERIAL
 cd $BUILD_SERIAL
 make ExperimentalUpdate
 make ExperimentalConfigure
 make ExperimentalBuild
 make ExperimentalSubmit

 ls -lt $DIRAC/$BUILD_SERIAL/dirac.x
 ldd $DIRAC/$BUILD_SERIAL/dirac.x
 python binary_info.py  > VERSION_cmake  2>&1  # get version

timestamp1=`date +\%F_\%k-\%M-\%S`; timestamp=${timestamp1// /};
echo -e "\n build 2 - Intel,MKL,serial,static - finished at $timestamp \n"

#################################################
### continue with the launching for the grid  ###
#################################################

# create the reduced tar-file containing all important stuff ...
cd $DIRAC

# copy there the mpirun static with related directories & files
cp    /home/ilias/bin/bin/openmpi_intel_stat/bin/mpirun                     $DIRAC/$BUILD/.
/bin/rm -rf $BUILD/share
mkdir $BUILD/share
cp -R -P /home/ilias/bin//home/ilias/bin/openmpi_intel_stat/share/openmpi   $DIRAC/$BUILD/share/.
cp -R  /home/ilias/bin/openmpi_intel_stat/etc                               $DIRAC/$BUILD/.

# problems with "--exclude test/*/*.out"  - so remove them as such
/bin/rm test/*/*.out*
/bin/rm -r test/*/compare

packed_dirac=dirac_grid_suite.tgz
echo -e "\n\n packing (slim) DIRAC suite for grid computations..."
tar czf $packed_dirac  pam  runtest VERSION test  basis  basis_dalton basis_ecp  $BUILD/dirac.x $BUILD/VERSION_cmake $BUILD/etc $BUILD/share $BUILD/mpirun  $BUILD_SERIAL/dirac.x $BUILD_SERIAL/VERSION_cmake  maintenance/grid_computing/diracrc_grid 
echo -e "\n packing done, we have tar-ball  file ready to be stored for grid-computing : "
ls -lt $packed_dirac
echo -e "\n Done. Tarball with static execs ready."

# exit here, don't do other stuff
echo -e "leaving here, you will have to upload tarball onto SE"
exit 0
