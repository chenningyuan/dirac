#!/bin/sh
#
# short script for copying and launching the dirac.x 
export rundir=`pwd`
echo "rundir=$rundir"
lcg-cp lfn://grid/compchem/ilias/dirac.x   file://$rundir/dirac.x  ||  { echo "Exiting...no dirac.x found..." ; exit 2; }
lcg-cp lfn://grid/compchem/ilias/pam       file://$rundir/pam      ||  { echo "Exiting...no pam found..."     ; exit 3; }
##########################################
##  make it executable and finally launch
##########################################
chmod 755 dirac.x
python ./pam --inp=fsccsd_IH.inp --mol=Mg.mol  --mw=92 --scratch=$PWD
