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
export PATH=/home/ilias/bin/openmpi-1.8.3_intel_static/bin:$PATH
export LD_LIBRARY_PATH=/home/ilias/bin/openmpi-1.8.3_intel_static/lib:$LD_LIBRARY_PATH
# my STATIC cmake stuff
export PATH=/home/ilias/bin/cmakestatic/bin:$PATH
echo "My own static CMake is in PATH."

# my own binutils-gold
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

cd $DIRAC

# submodule
git submodule init
git submodule update

export MATH_ROOT=$HOME/intel/mkl; echo "Activated MATH_ROOT=$MATH_ROOT"

# update the local trunk...
#echo -e "\n updating the local $DIRAC"
#git pull

BUILD_OMPI_INTEL=build_intelmkl_openmpi_i8_static
if [[ -d "$BUILD_OMPI_INTEL" ]]; then
   echo "removing previous build directory $BUILD_OMPI_INTEL"
  /bin/rm -rf $BUILD_OMPI_INTEL
fi
# deactivate interest module !
./pamadm -a --minus --interest --builddir=$BUILD_OMPI_INTEL
./setup --run-unit-tests  --fc=mpif90 --cc=mpicc --cxx=mpiCC --static -D BUILDNAME="grid_savba_ompi_intel_mkl_i8_STATIC" --int64 -D DART_TESTING_TIMEOUT=99999  $BUILD_OMPI_INTEL
 cd $BUILD_OMPI_INTEL
 ctest -D ExperimentalUpdate   --track Grid
 ctest -D ExperimentalConfigure   --track Grid
 ctest -D ExperimentalBuild  --track Grid
 ctest  --track Grid  -D ExperimentalTest -R cosci_energy 
 ctest -D ExperimentalSubmit  --track Grid

 ls -lt $DIRAC/$BUILD_OMPI_INTEL/dirac.x
 ldd $DIRAC/$BUILD_OMPI_INTEL/dirac.x
 python binary_info.py > VERSION_cmake  2>&1  # get version

timestamp1=`date +\%F_\%k-\%M-\%S`; timestamp=${timestamp1// /};
echo -e "\n build 1 - Intel-OpenMPI,MKL,i8,static - finished at $timestamp \n"

sleep 30

###################################################################################
cd $DIRAC
BUILD_SERIAL_INTEL=build_intelmkl_i8_static
if [[ -d "$BUILD_SERIAL_INTEL" ]]; then
   echo "removing previous build directory $BUILD_SERIAL_INTEL"
  /bin/rm -rf $BUILD_SERIAL_INTEL
fi

./pamadm -a --builddir=$BUILD_SERIAL_INTEL
./setup --fc=ifort --cc=icc --cxx=icpc --static -D BUILDNAME="grid_savba_intel_mkl_i8_STATIC" --int64 -D DART_TESTING_TIMEOUT=99999  $BUILD_SERIAL_INTEL
 cd $BUILD_SERIAL_INTEL
 ctest -D ExperimentalUpdate   --track Grid
 ctest -D ExperimentalConfigure   --track Grid
 ctest -D ExperimentalBuild    --track Grid
 ctest  --track Grid   -D ExperimentalTest -R cosci_energy  
 ctest -D ExperimentalSubmit  --track Grid

 ls -lt $DIRAC/$BUILD_SERIAL_INTEL/dirac.x
 ldd $DIRAC/$BUILD_SERIAL_INTEL/dirac.x
 python binary_info.py  > VERSION_cmake  2>&1  # get version

timestamp1=`date +\%F_\%k-\%M-\%S`; timestamp=${timestamp1// /};
echo -e "\n build 2 - Intel,MKL,i8,serial,static - finished at $timestamp \n"


sleep 30

###################################################################################
cd $DIRAC
BUILD_SERIAL_INTEL1=build_intel_ownmath_i8_static
if [[ -d "$BUILD_SERIAL_INTEL1" ]]; then
   echo "removing previous build directory $BUILD_SERIAL_INTEL1"
  /bin/rm -rf $BUILD_SERIAL_INTEL1
fi
./pamadm -a --builddir=$BUILD_SERIAL_INTEL1
./setup --fc=ifort --cc=icc --cxx=icpc  --blas=builtin  --lapack=builtin --static -D BUILDNAME="grid_savba_intel_ownmath_i8_STATIC" --int64 -D DART_TESTING_TIMEOUT=99999  $BUILD_SERIAL_INTEL1
 cd $BUILD_SERIAL_INTEL1
 ctest -D ExperimentalUpdate   --track Grid
 ctest -D ExperimentalConfigure   --track Grid
 ctest -D ExperimentalBuild    --track Grid
 ctest  --track Grid  -D ExperimentalTest -R cosci_energy  
 ctest -D ExperimentalSubmit  --track Grid

 ls -lt $DIRAC/$BUILD_SERIAL_INTEL1/dirac.x
 ldd $DIRAC/$BUILD_SERIAL_INTEL1/dirac.x
 python binary_info.py  > VERSION_cmake  2>&1  # get version
timestamp1=`date +\%F_\%k-\%M-\%S`; timestamp=${timestamp1// /};
echo -e "\n build 2.1 - Intel,ownmath,i8,serial,static - finished at $timestamp \n"

sleep 30

#####################
##  GNU+MKL STATIC ##
#####################
cd $DIRAC
BUILD_SERIAL_GNU=build_gnumkl_i8_static
if [[ -d "$BUILD_SERIAL_GNU" ]]; then
   echo "removing previous build directory $BUILD_SERIAL_GNU"
  /bin/rm -rf $BUILD_SERIAL_GNU
fi
./pamadm -a --builddir=$BUILD_SERIAL_GNU
./setup --fc=gfortran --cc=gcc --cxx=g++ --static -D BUILDNAME="grid_savba_gnu_mkl_i8_STATIC" --int64 -D DART_TESTING_TIMEOUT=99999  $BUILD_SERIAL_GNU
 cd $BUILD_SERIAL_GNU
 ctest -D ExperimentalUpdate      --track Grid
 ctest -D ExperimentalConfigure   --track Grid
 ctest -D ExperimentalBuild        --track Grid
 ctest --track Grid  -D ExperimentalTest -R cosci_energy  
 ctest -D ExperimentalSubmit  --track Grid

 ls -lt $DIRAC/$BUILD_SERIAL_GNU/dirac.x
 ldd $DIRAC/$BUILD_SERIAL_GNU/dirac.x
 python binary_info.py  > VERSION_cmake  2>&1  # get version

timestamp1=`date +\%F_\%k-\%M-\%S`; timestamp=${timestamp1// /};
echo -e "\n build 3 - GNU,MKL,i8,serial,static - finished at $timestamp \n"

sleep 30

#########################
##  GNU+OWNMATH STATIC ##
#########################
cd $DIRAC
BUILD_SERIAL_GNU1=build_gnu_ownmath_i8_static
if [[ -d "$BUILD_SERIAL_GNU1" ]]; then
   echo "removing previous build directory $BUILD_SERIAL_GNU1"
  /bin/rm -rf $BUILD_SERIAL_GNU1
fi
./pamadm -a --builddir=$BUILD_SERIAL_GNU1
./setup --fc=gfortran --cc=gcc --cxx=g++ --static --lapack=builtin --blas=builtin -D BUILDNAME="grid_savba_serial_gnu_ownmath_i8_STATIC" --int64 -D DART_TESTING_TIMEOUT=99999  $BUILD_SERIAL_GNU1
 cd $BUILD_SERIAL_GNU1
 ctest -D ExperimentalUpdate     --track Grid
 ctest -D ExperimentalConfigure   --track Grid
 ctest -D ExperimentalBuild   --track Grid
 ctest -D ExperimentalTest -R cosci_energy    --track Grid
 ctest -D ExperimentalSubmit  --track Grid    --track Grid

 ls -lt $DIRAC/$BUILD_SERIAL_GNU1/dirac.x
 ldd $DIRAC/$BUILD_SERIAL_GNU1/dirac.x
 python binary_info.py  > VERSION_cmake  2>&1  # get version

timestamp1=`date +\%F_\%k-\%M-\%S`; timestamp=${timestamp1// /};
echo -e "\n build 4 - GNU,ownmath,i8,serial,static - finished at $timestamp \n"

sleep 30

#################################################
### continue with the packing for grid servers ##
#################################################

# create the reduced tar-file containing all important stuff ...
cd $DIRAC

# copy there the mpirun static with related directories & files
cp    /home/ilias/bin/openmpi_intel_stat/bin/mpirun                         $DIRAC/$BUILD_OMPI_INTEL/.
/bin/rm -rf $BUILD_OMPI_INTEL/share
mkdir $BUILD_OMPI_INTEL/share
cp -R -P  /home/ilias/bin/openmpi_intel_stat/share/openmpi                  $DIRAC/$BUILD_OMPI_INTEL/share/.
cp -R     /home/ilias/bin/openmpi_intel_stat/etc                            $DIRAC/$BUILD_OMPI_INTEL/.
# don't forget export OPAL_PREFIX=$PWD/$BUILD 

# problems with "--exclude test/*/*.out"  - so remove them as such
/bin/rm    test/*/*.out
/bin/rm -r test/*/compare

packed_dirac=DIRAC_grid_suite.tgz
echo -e "\n\n packing (slim, with tests) DIRAC suite for grid computations..."
tar czf $packed_dirac  pam  VERSION test  basis  basis_dalton  basis_ecp \
$BUILD_OMPI_INTEL/dirac.x $BUILD_OMPI_INTEL/VERSION_cmake $BUILD_OMPI_INTEL/etc $BUILD_OMPI_INTEL/share $BUILD_OMPI_INTEL/mpirun  \
$BUILD_SERIAL_INTEL/dirac.x $BUILD_SERIAL_INTEL/VERSION_cmake \
$BUILD_SERIAL_INTEL1/dirac.x $BUILD_SERIAL_INTEL1/VERSION_cmake \
$BUILD_SERIAL_GNU/dirac.x $BUILD_SERIAL_GNU/VERSION_cmake \
$BUILD_SERIAL_GNU1/dirac.x $BUILD_SERIAL_GNU1/VERSION_cmake \
/home/ilias/bin/cmakestatic \
maintenance/grid_computing/diracrc_grid 

echo -e "\n packing done, we have tar-ball file ready to be placed on SE for grid-computing : "
ls -lt $packed_dirac
echo -e "\n Done. Tarball with static execs ready."

# exit here, don't do other stuff
echo -e "\n Leaving here, you will have to upload the tarball manually onto your SE."
exit 0
