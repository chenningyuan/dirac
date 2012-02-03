#!/bin/sh
#
# short script for copying and launching the DIRAC testing machinery

echo -e "\n --- Hostname: \c"; hostname -f; echo
echo -e "I am here: $PWD \n"
echo -e "My space: \c"; df -h .

echo "files on my voce space lfn://grid/voce/ilias:"
lcg-ls lfn://grid/voce/ilias

###  copy files to the current CE ###
export rundir=`pwd`
lcg-cp lfn://grid/voce/ilias/dirac_test_suite.tgz  file://$rundir/dirac_test_suite.tgz  ||  { echo "Exiting...no dirac_test_suite.tgz found..." ; exit 3; }

##########################################
##  make it executable and finally launch
##########################################
echo -e "\n Job started at : \c"; date

 export OPAL_PREFIX=$PWD/build
 echo -e "OPAL_PREFIX=$OPAL_PREFIX"
 
 export PATH=$PWD/build:$PATH
 echo -e "PATH=$PATH"
 echo -e "LD_LIBRARY_PATH=$LD_LIBRARY_PATH"

echo -e "\n ---- server environment ----"
set

echo -e "\n ---- server CPU equipment ----"
cat /proc/cpuinfo

echo -e "\n ---- server memory equipment ----"
cat /proc/meminfo

echo -e "\n ---- DIRAC unpacking ---- \n"
tar -xvf dirac_test_suite.tgz
chmod  u+x build/dirac.x build/mpirun pam runtest

# restore diracrc file for paralle run on the grid
cp diracrc_grid diracrc

echo -e "\n --- Launching simple pam test  --- \n "; 

python ./pam --inp=test/fscc/fsccsd_IH.inp --mol=test/fscc/Mg.mol  --mw=92 --scratch=$PWD --outcmo --mpi=2

#####################################################################

echo -e "\n --- Launching runtest --- \n "; 

#python ./runtest --verbose --tests=quick --mpi=2   # works ...
python ./runtest --verbose --tests="cosci_energy fscc acmoin" --mpi=2 # works
#python ./runtest --verbose --quick # works
#python ./runtest --verbose --all  # works

############################
#  print out some info .. 
############################
echo -e "\n -------------------------------------------------------"
echo -e "Job finished at: \c"; date
echo -e "My space: \c"; df -h .
echo -e "I am here: \c";pwd
echo -e  "All files:\n"; ls -lt
echo -e "Total occupied space: \c"; du -h -s .
