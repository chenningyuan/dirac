#!/bin/sh
#
# short script for copying and launching the DIRAC testing machinery

export rundir=`pwd`

echo -e "Hostname: \c"; hostname -f
echo -e "I am here: \c";$rundir

echo "files on my voce space lfn://grid/voce/ilias:"
lcg-ls lfn://grid/voce/ilias

###  copy files to the current CE ###
#lcg-cp lfn://grid/compchem/ilias/dirac.x   file://$rundir/dirac.x  ||  { echo "Exiting...no dirac.x found..." ; exit 2; }
#lcg-cp lfn://grid/compchem/ilias/pam       file://$rundir/pam      ||  { echo "Exiting...no pam found..."     ; exit 3; }

lcg-cp lfn://grid/voce/ilias/dirac_test_suite.tgz    file://$rundir/dirac_test_suite.tgz  ||  { echo "Exiting...no dirac_test_suite.tgz found..." ; exit 3; }

##########################################
##  make it executable and finally launch
##########################################
echo -e "\n Job started at : \c"; date

tar -xvf dirac_test_suite.tgz
chmod  u+x build/dirac.x pam runtest

echo -e "\n --- Launching runtest --- \n "; 

#python ./pam --inp=fsccsd_IH.inp --mol=Mg.mol  --mw=92 --scratch=$PWD --outcmo
#python ./runtest --verbose --tests=cosci_energy # works
#python ./runtest --verbose --tests="cosci_energy fscc acmoin" # works
#python ./runtest --verbose --quick # works
python ./runtest --verbose --all

############################
#  print out some info .. 
############################
echo "\n -------------------------------------------------------"
echo -e "Job finished at: \c"; date
echo -e "My space: \c"; df -h .
echo -e "I am here: \c";pwd
echo  "All files:"; ls -lt
echo -e "Total occupied space: \c"; du -h -s .
