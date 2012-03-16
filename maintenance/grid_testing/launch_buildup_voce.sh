#!/bin/sh
#######################################################################################################
#
#
#               short script for copying and launching the DIRAC testing machinery
#
#
#######################################################################################################


echo -e "\n --- Hostname: \c"; hostname -f; echo
echo -e " ---------  I am here: --------- \n $PWD "
echo -e "\n---  My space: "; df -h .

echo -e  "\n --- files on my vo=voce space lfn://grid/voce/ilias:"
lcg-ls -l lfn://grid/voce/ilias

#####################################
###  copy files to the current CE ###
#####################################
export rundir=`pwd`
lcg-cp lfn://grid/voce/ilias/dirac_current.tgz  file://$rundir/dirac_current.tgz 

##########################################
##  make it executable and finally launch
##########################################
echo -e "\n Job started at : \c"; date

export OPAL_PREFIX=$PWD/static_execs
echo -e "OPAL_PREFIX=$OPAL_PREFIX"
 
export PATH=$OPAL_PREFIX:$PATH
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
tar xvfz dirac_current.tgz
chmod u+x static_execs/mpirun static_execs/cmake static_execs/ctest  pam  runtest setup

echo -e "\n ---  Launching pamadm configuration  --- \n "; 
python pamadm -a

echo -e "\n --- Launching setup configuration  --- \n"; 
python  ./setup  -D VERBOSE_OUTPUT=ON --auto --debug

echo -e "\n --- Launching CMake buildup configuration  --- \n"; 
cd build
cmake ..
make dirac.x

echo -e "\n --- Launching simple pam test(s)  --- \n "; 
cd ...
#python ./pam --inp=test/fscc/fsccsd_IH.inp --mol=test/fscc/Mg.mol  --mw=92 --scratch=$PWD --outcmo --mpi=2
python ./pam --inp=test/fscc/fsccsd_IH.inp --mol=test/fscc/Mg.mol  --mw=92  

#####################################################################
echo -e "\n --- Launching runtest --- \n "; 
#python ./runtest --verbose --tests=quick --mpi=2   # works ...
#python ./runtest --verbose --tests="cosci_energy fscc acmoin" --mpi=2 # works
python ./runtest --verbose --tests="cosci_energy" 
#python ./runtest --verbose --quick --mpi=2  # works
#python ./runtest --verbose --all  # works

#############################################################
#
# own cmake stuff
# cd build
# cmake  -DVERBOSE_OUTPUT=ON -DBUILDNAME=grid_test -DENABLE_MPI=ON -DRUNTEST_NR_PROC=2 ..
# make ExperimentalStart
# make ExperimentalTest
# make ExperimentalSubmit
#

echo -e "\n --- Want packing all wanted stuff from the grid CE  --- \n "; 
ls -lt
tar --version
#tar --help
#echo " "
#tar cvzf dirac_test_suite_back.tgz test *.out *.tgz

############################
#  print out some info .. 
############################
echo -e "\n -------------------------------------------------------"
echo -e "Job finished at: \c"; date
echo -e "I am here: \c"; pwd
echo -e "My space: "; df -h .
echo -e "\n All files:"; ls -lt
echo -e "Total occupied space: \c"; du -h -s .
