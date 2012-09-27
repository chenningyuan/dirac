#!/bin/sh
#######################################################################################################
#
#
#           Script for copying and launching the DIRAC suite on grid CE
#
#
#######################################################################################################

# check the parameter - VO name

if [[ $1 != "voce" && $1 != "compchem" && $1 != "isragrid" && $1 != "osg" ]]; then
 echo -e "\n wrong parameter - VO name : $1 "
 exit -12
else
 VO=$1
 echo -e "OK, you specified properly the VO=$VO"
 #exit 0
fi

THIS_HOSTNAME=`hostname -f`
echo -e "\n --- Hostname: $THIS_HOSTNAME "; 
echo -e " ---------  My current directory :  $PWD "
echo -e "---  My space, df -h: "; df -h .

echo -e "-- \n Is there cmake ? which cmake="; which cmake

echo -e  "\n Files on my vo=$VO space, lcg-ls -l lfn://grid/$VO/ilias:"
lcg-ls -l lfn://grid/$VO/ilias
echo -e "\n  File attributes lcg-lr  lfn://grid/$VO/ilias/dirac_grid_suite.tgz:"
lcg-lr    lfn://grid/$VO/ilias/dirac_grid_suite.tgz
echo " ------------------------------------------------------ "

#####################################
###  download file to the current CE 
#####################################
echo -e "\n Going to download from SE, lcg-cp lfn://grid/$VO/ilias/dirac_grid_suite.tgz  file://$PWD/dirac_grid_suite.tgz...."
lcg-cp -v lfn://grid/$VO/ilias/dirac_grid_suite.tgz   file://$PWD/dirac_grid_suite.tgz || { echo "Error exit,downloading of dirac_grid_suite.tg failed" ; exit -2; }
echo -e "\nDownloaded ? File attributes (ls -lt dirac_grid_suite.tgz) :"
ls -lt dirac_grid_suite.tgz || { echo -e "No dirac_grid_suite.tgz in local directory found, error exit"; exit -3;}
echo " ---------------------------------------------------------------------------------------------------------------------"

##########################################
##  make it executable and finally launch
##########################################
echo -e "\n Job started at : \c"; date

echo -e "Server specification: \c";  uname -a

echo -e "\n server Python version: \c"; which python
python -V
echo -e "\n"

# directory with all static executables - dirac.x, cmake, ctest, mpirun - and related directories
 export BUILD=build_Intel_OpenMPI_MKL_ILP64_STATIC
 export BUILD_SERIAL=build_Intel_serial_MKL_ILP64_STATIC

 unset OPAL_PREFIX
 export OPAL_PREFIX=$PWD/$BUILD
 echo -e "Variable OPAL_PREFIX=$OPAL_PREFIX"
 export PATH=$PWD/$BUILD:$PATH
 #export PATH=$PWD/$BUILD_SERIAL:$PWD/$BUILD:$PATH
 echo -e "Modified PATH=$PATH"
 echo -e "LD_LIBRARY_PATH=$LD_LIBRARY_PATH"

#echo -e "\n ---- server environment ---- set :"
#set

#echo -e "\n ---- server CPU equipment, /proc/cpuinfo ----"
#cat /proc/cpuinfo

# get number of CPUs on the given node
nprocs=`cat /proc/cpuinfo | grep processor | wc -l`
#echo -e "\n Number of CPUs :\c";cat /proc/cpuinfo | grep processor | wc -l
echo -e "\n Number of CPUs : $nprocs";
echo -e " CPUs model : \c"; cat /proc/cpuinfo | grep "model name"  | tail -n 1

echo -e "\n ---- server total memory equipment ----"
cat /proc/meminfo | grep MemTotal

echo -e "\n ---- DIRAC unpacking with tar : \c "; which tar;  tar --version
tar xfz dirac_grid_suite.tgz

echo -e "\n Content of directory BUILD=$BUILD: \c"; ls -lt $BUILD
echo -e "\n Content of directory BUILD_SERIAL=$BUILD_SERIAL: \c"; ls -lt $BUILD_SERIAL

#echo -e "\n making executables:"
#ls -lt $BUILD/dirac.x  $BUILD/cmake  $BUILD/ctest $BUILD/mpirun  pam  runtest pamadm setup
#chmod  u+x  $BUILD/dirac.x  $BUILD/cmake  $BUILD/ctest $BUILD/mpirun  pam  runtest pamadm setup
#ls -lt $BUILD/dirac.x  $BUILD/cmake  $BUILD/ctest $BUILD/mpirun  pam  runtest pamadm setup

# create diracrc file for paralle run on the grid CE
cp maintenance/grid_computing/diracrc_grid  diracrc
echo -e "\n ---- diracrc created: "; cat diracrc

echo -e "\n --- Launching simple pam test  --- \n "; 
#python ./pam --inp=test/fscc/fsccsd_IH.inp --mol=test/fscc/Mg.mol  --mw=92 --outcmo --mpi=2
python ./pam --inp=test/fscc/fsccsd_IH.inp --mol=test/fscc/Mg.mol  --mw=92 --outcmo --mpi=$nprocs
#python ./pam --inp=test/fscc/fsccsd_IH.inp --mol=test/fscc/Mg.mol  --mw=92 --outcmo 
#python ./pam --inp=test/fscc/fsccsd_IH.inp --mol=test/fscc/Mg.mol  --mw=92 --scratch=$PWD --outcmo 

#####################################################################

echo -e "\n --- Launching quick runtest --- \n "; 

#python ./runtest --verbose --tests=quick --mpi=2   # works ...
#python ./runtest --verbose --tests="cosci_energy fscc acmoin" --mpi=2 # works
python ./runtest --verbose --quick --mpi=$nprocs  #  works
#python ./runtest --verbose --quick #  works
#python ./runtest --verbose --all  # works

#############################################################
#               check own static cmake stuff
#############################################################

 python ./pamadm -a

 cp maintenance/grid_computing/diracrc_build_grid  $BUILD/diracrc
 cp maintenance/grid_computing/diracrc_build_grid  $BUILD_SERIAL/diracrc
 echo -e "\n ---- own diracrc in build-directories created: "; cat $BUILD/diracrc

 #cd $BUILD_SERIAL
 cd $BUILD
 echo -e "\n I am in build-directory  $PWD \n"
 echo -e "--- Directory $BUILD content ---:"; ls -lt

echo -e "\n False configuration in order to enable tests running"
# FIX ERROR:
# CMake Error: The current CMakeCache.txt directory /tmp/ilias_OUWzTUf41_ntFwYFY3Pk7Q/build/CMakeCache.txt is different than the directory /home/ilias/qch_work/software/dirac_git/trunk_current/build where CMakeCache.txt was created. This may result in binaries being created in the wrong place. If you are not sure, reedit the CMakeCache.txt
#CMake Error: The source "/tmp/ilias_OUWzTUf41_ntFwYFY3Pk7Q/CMakeLists.txt" does not match the source "/home/ilias/qch_work/software/dirac_git/trunk_current/CMakeLists.txt" used to generate cache.  Re-run cmake with a different source directory.
mv CMakeCache.txt CMakeCache.txt_STORED

#note: make  -DENABLE_MPI=OFF in order to avoid cmake's MPI configuration on an unknown grid's CE
 ./cmake  -DVERBOSE_OUTPUT=ON -DBUILDNAME="$THIS_HOSTNAME.grid_test_ompi_intel_mkl_static" -DENABLE_MPI=OFF  -DRUNTEST_NR_PROC=$nprocs -DRUNTEST_TIMEOUT="00:18:00" -DDART_TESTING_TIMEOUT=99999 -D RUNTEST_BINARY=$DIRAC/$BUILD/dirac.x ..

# make ExperimentalConfigure
 make ExperimentalStart
 make ExperimentalTest
 make ExperimentalSubmit

## return back into main directory
 cd ..
 echo -e "\n I am in main directory $PWD \n"

# ----------------------------------------- #
# ---- try own buildup on the grid CE ----- #
# ----------------------------------------- #
 echo -e "\n pgf90 ? \c"; which pgf90
 echo -e "\n ifort ? \c"; which ifort
 echo -e "\n echo MKLROOT=$MKLROOT"
 echo -e "\n lfc ? \c"; which lfc
 echo -e "\n gfortran44 ? \c"; which gfortran44
 echo -e "\n gfortran ? \c"; which gfortran
 echo -e "\n g++ ? \c"; which g++
 echo -e "\n mpif90 ? \c"; which mpif90
 echo -e "\n mpicc ? \c"; which mpicc
 echo -e "\n mpiCC ? \c"; which mpiCC

 echo -e "\n ==== Try own buildup of the host server in build directory ==== "

 #python ./setup -D VERBOSE_OUTPUT=ON --fc=$FORTRAN
 python ./setup -D VERBOSE_OUTPUT=ON --int64=on --mpi=off -D BUILDNAME="$THIS_HOSTNAME_grid_ownbuild" -D RUNTEST_TIMEOUT="00:40:00" -D DART_TESTING_TIMEOUT=99999
 cd build
 cmake ..
 make Experimental


# --------------------------------------------------------------------
echo -e "\n --------------------------------- \n "; 
# return to the "main" directory...
cd  ..
# delete old tar-ball first
ls -lt dirac_grid_suite.tgz
echo -e "\n deleting the old dirac_grid_suite.tgz..."
rm dirac_grid_suite.tgz
echo "check files..."
ls -lt
echo -e "\n --- Packing all wanted stuff from the grid CE  --- "; 
tar --version
#tar --help
echo " "
#tar cvzf dirac_suite_back.tgz test *.out *.tgz
#tar cvzf dirac_suite_back.tgz test *.out 
tar czf dirac_grid_suite_back.tgz  *

############################
#  print out some info .. 
############################
echo -e "\n -------------------------------------------------------"
echo -e "Job finished at: \c"; date
echo -e "I am here: \c"; pwd
echo -e "My space: "; df -h .
echo -e "\n All files:"; ls -lt
echo -e "\n Total occupied space: \c"; du -h -s .

exit 0
