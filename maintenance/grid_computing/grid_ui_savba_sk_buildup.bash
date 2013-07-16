#!/bin/bash
#############################################################################
#
#      parametrized script for DIRAC buildup in the grid environment
#
#############################################################################

##################################################
# fill the PATH to be the same as in online mode
##################################################
source $HOME/intel/bin/compilervars.sh intel64
echo "Intel Fortran/C/C++ noncommercial compilers with MKL library activated."
# My own OpenMPI Intel static
export PATH=/home/ilias/bin/openmpi_intel_stat/bin:$PATH
export LD_LIBRARY_PATH=/home/ilias/bin/openmpi_intel_stat/lib:$LD_LIBRARY_PATH
# my cmake stuff
export PATH=/home/ilias/bin/newer_cmake:$PATH
echo "My own (dynamic) CMake is in PATH."

# my binutils-gold
export LIBDIR=/home/ilias/bin/binutils-gold/lib:$LIBDIR
export LD_LIBRARY_PATH=$LIBDIR:$LD_LIBRARY_PATH
export LD_RUN_PATH=$LIBDIR:$LD_RUN_PATH
export PATH=/home/ilias/bin/binutils-gold/binutils-master/gold:$PATH

echo "My PATH=$PATH"
#echo "which make:"; which make
echo -e "\n set:"; set

timestamp1=`date +\%F_\%k-\%M-\%S`; timestamp=${timestamp1// /};
echo "Running DIRAC cdash buildup at "$timestamp

DIRAC=/home/ilias/Work/qch/software/trunk
#
BUILD_OMPI_INTEL=build_intel_openmpi_static
BUILD_SERIAL_INTEL=build_serial_intel_static
BUILD_SERIAL_GNU=build_serial_gnu_static

cd $DIRAC

export MATH_ROOT=$HOME/intel/mkl; echo "Activated MATH_ROOT=$MATH_ROOT"

# update the local trunk...
echo -e "\n updating the local $DIRAC"
git pull

if [[ -d "$BUILD_OMPI_INTEL" ]]; then
   echo "removing previous build directory $BUILD_OMPI_INTEL"
  /bin/rm -rf $BUILD_OMPI_INTEL
fi
./setup --fc=mpif90 --cc=mpicc --cxx=mpiCC --static -D BUILDNAME="grid_savba_ompi_intel_mkl_ilp64_STATIC" -D RUNTEST_NR_PROC=2 --int64 -D RUNTEST_TIMEOUT="00:18:00" -D DART_TESTING_TIMEOUT=99999  $BUILD_OMPI_INTEL
 cd $BUILD_OMPI_INTEL
 make ExperimentalUpdate
 make ExperimentalConfigure
 make ExperimentalBuild
 make ExperimentalSubmit

 ls -lt $DIRAC/$BUILD_OMPI_INTEL/dirac.x
 ldd $DIRAC/$BUILD_OMPI_INTEL/dirac.x
 python binary_info.py > VERSION_cmake  2>&1  # get version

timestamp1=`date +\%F_\%k-\%M-\%S`; timestamp=${timestamp1// /};
echo -e "\n build 1 - Intel-OpenMPI,MKL,static - finished at $timestamp \n"

cd $DIRAC
if [[ -d "$BUILD_SERIAL_INTEL" ]]; then
   echo "removing previous build directory $BUILD_SERIAL_INTEL"
  /bin/rm -rf $BUILD_SERIAL_INTEL
fi

./setup --fc=ifort --cc=icc --cxx=icpc --static -D BUILDNAME="grid_savba_serial_intel_mkl_ilp64_STATIC" --int64 -D RUNTEST_TIMEOUT="00:18:00" -D DART_TESTING_TIMEOUT=99999  $BUILD_SERIAL_INTEL
 cd $BUILD_SERIAL_INTEL
 make ExperimentalUpdate
 make ExperimentalConfigure
 make ExperimentalBuild
 make ExperimentalSubmit

 ls -lt $DIRAC/$BUILD_SERIAL_INTEL/dirac.x
 ldd $DIRAC/$BUILD_SERIAL_INTEL/dirac.x
 python binary_info.py  > VERSION_cmake  2>&1  # get version

timestamp1=`date +\%F_\%k-\%M-\%S`; timestamp=${timestamp1// /};
echo -e "\n build 2 - Intel,MKL,serial,static - finished at $timestamp \n"

#################################################
### continue with the launching for the grid  ###
#################################################

# create the reduced tar-file containing all important stuff ...
cd $DIRAC

# copy there the mpirun static with related directories & files
cp    /home/ilias/bin/bin/openmpi_intel_stat/bin/mpirun                     $DIRAC/$BUILD_OMPI_INTEL/.
/bin/rm -rf $BUILD_OMPI_INTEL/share
mkdir $BUILD_OMPI_INTEL/share
cp -R -P /home/ilias/bin//home/ilias/bin/openmpi_intel_stat/share/openmpi   $DIRAC/$BUILD_OMPI_INTEL/share/.
cp -R  /home/ilias/bin/openmpi_intel_stat/etc                               $DIRAC/$BUILD_OMPI_INTEL/.

# problems with "--exclude test/*/*.out"  - so remove them as such
/bin/rm test/*/*.out*
/bin/rm -r test/*/compare

packed_dirac=DIRAC_grid_suite.tgz
echo -e "\n\n packing (slim) DIRAC suite for grid computations..."
tar czf $packed_dirac  pam  runtest VERSION test  basis  basis_dalton basis_ecp  $BUILD_OMPI_INTEL/dirac.x $BUILD_OMPI_INTEL/VERSION_cmake $BUILD_OMPI_INTEL/etc $BUILD_OMPI_INTEL/share $BUILD_OMPI_INTEL/mpirun  $BUILD_SERIAL_INTEL/dirac.x $BUILD_SERIAL_INTEL/VERSION_cmake  maintenance/grid_computing/diracrc_grid 
echo -e "\n packing done, we have tar-ball  file ready to be placed on SE for grid-computing : "
ls -lt $packed_dirac
echo -e "\n Done. Tarball with static execs ready."

# exit here, don't do other stuff
echo -e "leaving here, you will have to upload tarball onto SE"
exit 0
