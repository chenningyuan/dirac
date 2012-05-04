#!/bin/sh
#
# short script for copying and launching the dirac.x 
#

# first copy executable to the grid storage space
DIRAC=/home/ilias/qch_work/software/dirac_git/trunk/build

### executable must be present in the current directory!
rundir=`pwd`
lcg-cp lfn://grid/voce/ilias/dirac.x file://$rundir/dirac.x  || { echo "Exiting...no dirac.x found..." ; exit 2; }
#lcg-cp lfn://grid/voce/ilias/dirac.x file://$DIRAC/dirac.x  || { echo "Exiting...no dirac.x found..." ; exit 2; } 

# make it executable and finally launch
chmod 755 dirac.x
export DIRWRK=84000000
./dirac.x
