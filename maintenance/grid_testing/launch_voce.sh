#!/bin/sh
#
# short script for copying and launching the DIRAC testing machinery
export rundir=`pwd`
echo "rundir=$rundir"

echo "files on my voce space lfn://grid/voce/ilias:"
lcg-ls lfn://grid/voce/ilias

###  copy files to the current CE ###
#lcg-cp lfn://grid/compchem/ilias/dirac.x   file://$rundir/dirac.x  ||  { echo "Exiting...no dirac.x found..." ; exit 2; }
#lcg-cp lfn://grid/compchem/ilias/pam       file://$rundir/pam      ||  { echo "Exiting...no pam found..."     ; exit 3; }

lcg-cp lfn://grid/voce/ilias/dirac_test_suite.tgz    file://$rundir/dirac_test_suite.tgz  ||  { echo "Exiting...no dirac_test_suite.tgz found..." ; exit 3; }

##########################################
##  make it executable and finally launch
##########################################
tar -xvf dirac_test_suite.tgz
chmod  u+x build/dirac.x pam runtest

#python ./pam --inp=fsccsd_IH.inp --mol=Mg.mol  --mw=92 --scratch=$PWD --outcmo
#python ./runtest --verbose --tests=cosci_energy
python ./runtest --verbose --tests="cosci_energy fscc acmoin"

############################
#  print out some info .. 
############################
echo "Done, date, hostname, current directory, list of files:"
date
echo "Hostname:"; hostname -f
echo "My space:"; df -h .
pwd
ls -lt
echo "Total occupied space:"
du -h -s .
