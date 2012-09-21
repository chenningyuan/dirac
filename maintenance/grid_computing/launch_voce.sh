#!/bin/sh
#######################################################################################################
#
#
#           Script for copying and launching the DIRAC suite on grid
#
#
#######################################################################################################

echo -e "\n --- Hostname: \c"; hostname -f; echo
echo -e " ---------  My current directory :  $PWD "
echo -e "---  My space: "; df -h .

echo -e  "\n Files on my vo=voce space lfn://grid/voce/ilias:"
lcg-ls -l lfn://grid/voce/ilias
lcg-lr    lfn://grid/voce/ilias/dirac_grid_suite.tgz
echo " ------------------------- "

#####################################
###  download file to the current CE ###
#####################################
echo -e "\n Going to download lcg-cp lfn://grid/voce/ilias/dirac_grid_suite.tgz  file://$PWD/dirac_grid_suite.tgz...."
lcg-cp -v lfn://grid/voce/ilias/dirac_grid_suite.tgz  file://$PWD/dirac_grid_suite.tgz  
echo -e "\nDownloaded ? File attributes;  ls -lt dirac_grid_suite.tgz"
ls -lt dirac_grid_suite.tgz
echo " -----------------------"


##########################################
##  make it executable and finally launch
##########################################
echo -e "\n Job started at : \c"; date

echo -e "Server specification: \c";  uname -a

echo -e "\n server Python version: \c"; which python; python -V

echo -e "\n"

 # directory with all static executables - dirac.x, cmake, ctest, mpirun - and related directories
 export BUILD=build_Intel_OpenMPI_MKL_ILP64_STATIC
 export BUILD_SERIAL=build_Intel_serial_MKL_ILP64_STATIC

 export OPAL_PREFIX=$PWD/$BUILD
 echo -e "Variable OPAL_PREFIX=$OPAL_PREFIX"
 
 #export PATH=$PWD/$BUILD:$PATH
 export PATH=$PWD/$BUILD_SERIAL:$PWD/$BUILD:$PATH
 echo -e "Modified PATH=$PATH"
 echo -e "LD_LIBRARY_PATH=$LD_LIBRARY_PATH"

echo -e "\n ---- server environment ---- set :"
set

#echo -e "\n ---- server CPU equipment, /proc/cpuinfo ----"
#cat /proc/cpuinfo

echo -e "\n Number of CPUs :\c";cat /proc/cpuinfo | grep processor | wc -l

echo -e "\n ---- server memory equipment ----"
cat /proc/meminfo

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
python ./pam --inp=test/fscc/fsccsd_IH.inp --mol=test/fscc/Mg.mol  --mw=92 --outcmo 
#python ./pam --inp=test/fscc/fsccsd_IH.inp --mol=test/fscc/Mg.mol  --mw=92 --scratch=$PWD --outcmo 

#####################################################################

echo -e "\n --- Launching quick runtest --- \n "; 

#python ./runtest --verbose --tests=quick --mpi=2   # works ...
#python ./runtest --verbose --tests="cosci_energy fscc acmoin" --mpi=2 # works
#python ./runtest --verbose --quick --mpi=2  #  works
python ./runtest --verbose --quick #  works
#python ./runtest --verbose --all  # works

#############################################################
#               check own static cmake stuff
#############################################################

 python ./pamadm -a

 cd $BUILD_SERIAL
 echo -e "\n I am in $PWD \n"

# FIX ERROR:
# CMake Error: The current CMakeCache.txt directory /tmp/ilias_OUWzTUf41_ntFwYFY3Pk7Q/build/CMakeCache.txt is different than the directory /home/ilias/qch_work/software/dirac_git/trunk_current/build where CMakeCache.txt was created. This may result in binaries being created in the wrong place. If you are not sure, reedit the CMakeCache.txt
#CMake Error: The source "/tmp/ilias_OUWzTUf41_ntFwYFY3Pk7Q/CMakeLists.txt" does not match the source "/home/ilias/qch_work/software/dirac_git/trunk_current/CMakeLists.txt" used to generate cache.  Re-run cmake with a different source directory.
mv CMakeCache.txt CMakeCache.txt_STORED

# ./cmake  -DVERBOSE_OUTPUT=ON -DBUILDNAME="$HOSTNAME.grid_test_ompi_intel_mkl_static" -DENABLE_MPI=ON -DRUNTEST_NR_PROC=2 -DRUNTEST_TIMEOUT="00:18:00" -DDART_TESTING_TIMEOUT=99999 -DENABLE_SANITY_CHECKS=OFF -D RUNTEST_BINARY=$DIRAC/$BUILD/dirac.x ..
 ./cmake  -DVERBOSE_OUTPUT=ON -DBUILDNAME="$HOSTNAME.grid_serial_intel_mkl_static" -DENABLE_MPI=OFF  -DRUNTEST_TIMEOUT="00:18:00" -DDART_TESTING_TIMEOUT=99999 -DENABLE_SANITY_CHECKS=OFF -DRUNTEST_BINARY=$DIRAC/$BUILD_SERIAL/dirac.x ..

 make ExperimentalConfigure
 make ExperimentalStart
 make ExperimentalTest
 make ExperimentalSubmit
## return back into main directory
 cd ..

# ----------------------------------------- #
# ---- try own buildup on the grid CE ----- #
# ----------------------------------------- #
 echo -e "\n pgf90 ? \c"; which pgf90
 echo -e "\n ifort ? \c"; which ifort
 echo -e "\n gfortran44 ? \c"; which gfortran44
 echo -e "\n gfortran ? \c"; which gfortran

 echo -e "\n Own buildup of the host server "

 #python ./setup -D VERBOSE_OUTPUT=ON --fc=$FORTRAN
 python ./setup -D VERBOSE_OUTPUT=ON --int64=on --mpi=off
 cd build
 cmake ..
 make Experimental


# --------------------------------------------------------------------
echo -e "\n --- Packing all wanted stuff from the grid CE  --- \n "; 
# delete old tar-ball first
rm dirac_grid_suite.tgz
ls -lt
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
echo -e "Total occupied space: \c"; du -h -s .
