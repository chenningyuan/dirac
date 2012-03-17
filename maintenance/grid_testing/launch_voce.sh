#!/bin/sh
#######################################################################################################
#
#
#               short script for copying and launching the DIRAC testing machinery
#
#
#######################################################################################################

echo -e "\n --- Hostname: \c"; hostname -f; echo
echo -e " ---------  I am here: ---------  $PWD "
echo -e "---  My space: "; df -h .

echo -e  "\n files on my voce space lfn://grid/voce/ilias:"
lcg-ls -l lfn://grid/voce/ilias

echo " "

#####################################
###  copy files to the current CE ###
#####################################
export rundir=`pwd`
#lcg-cp lfn://grid/voce/ilias/dirac_test_suite.tgz  file://$rundir/dirac_test_suite.tgz  ||  { echo "Exiting...no dirac_test_suite.tgz found..." ; exit 3; }
lcg-cp lfn://grid/voce/ilias/dirac_current.tgz  file://$rundir/dirac_current.tgz  ||  { echo "Exiting...no dirac_current.tgz  found..." ; exit -1; }

##########################################
##  make it executable and finally launch
##########################################
echo -e "\n Job started at : \c"; date

 # directory with all static executables - dirac.x, cmake, ctest, mpirun - and related directories
 export BUILD=build_imported_static

 export OPAL_PREFIX=$PWD/$BUILD
 echo -e "OPAL_PREFIX=$OPAL_PREFIX"
 
 export PATH=$PWD/$BUILD:$PATH
 echo -e "PATH=$PATH"
 echo -e "LD_LIBRARY_PATH=$LD_LIBRARY_PATH"

echo -e "\n ---- server environment ----"
set

echo -e "\n ---- server CPU equipment ----"
cat /proc/cpuinfo

echo -e "\n Number of CPUs :\c";cat /proc/cpuinfo | grep processor | wc -l

echo -e "\n ---- server memory equipment ----"
cat /proc/meminfo

echo -e "\n ---- DIRAC unpacking ---- \n"
tar --version
tar xvfz dirac_test_suite.tgz
chmod  u+x $BUILD/dirac.x  $BUILD/mpirun $BUILD/cmake $BUILD/ctest pam runtest

# create diracrc file for paralle run on the grid
echo -e "\n ---- diracrc created ---- \n"
cp diracrc_grid diracrc

echo -e "\n --- Launching simple pam test  --- \n "; 

python ./pam --inp=test/fscc/fsccsd_IH.inp --mol=test/fscc/Mg.mol  --mw=92 --scratch=$PWD --outcmo --mpi=2

#####################################################################

echo -e "\n --- Launching runtest --- \n "; 

#python ./runtest --verbose --tests=quick --mpi=2   # works ...
#python ./runtest --verbose --tests="cosci_energy fscc acmoin" --mpi=2 # works
python ./runtest --verbose --quick --mpi=2  # works
#python ./runtest --verbose --all  # works

#############################################################
# own cmake stuff
#############################################################
 cd $BUILD
 cmake  -DVERBOSE_OUTPUT=ON -DBUILDNAME=grid_test -DENABLE_MPI=ON -DRUNTEST_NR_PROC=2 ..
 make ExperimentalStart
 make ExperimentalTest
 make ExperimentalSubmit

# ------------
echo -e "\n --- Packing all wanted stuff from the grid CE  --- \n "; 
ls -lt
tar --version
#tar --help
echo " "
tar cvzf dirac_test_suite_back.tgz test *.out *.tgz

############################
#  print out some info .. 
############################
echo -e "\n -------------------------------------------------------"
echo -e "Job finished at: \c"; date
echo -e "I am here: \c"; pwd
echo -e "My space: "; df -h .
echo -e "\n All files:"; ls -lt
echo -e "Total occupied space: \c"; du -h -s .
