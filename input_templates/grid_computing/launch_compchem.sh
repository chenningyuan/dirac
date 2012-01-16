#!/bin/sh
#
# short script for copying and launching the dirac.x 
#
#lcg-cr -d grid2.fe.infn.it  --vo compchem  -l lfn:/grid/compchem/ilias/dirac.x   dirac.x
#guid:2ef763fd-6265-40d3-b212-93e4f465476c
#ilias@grid3.ui.savba.sk:~/qch_work/software/dirac_git/trunk/input_templates/grid_computing/.lfc-ls /grid/compchem/ilias
#dirac.x
#dirac_space

# first copy executable to the grid storage space
#DIRAC=/home/ilias/qch_work/software/dirac_git/trunk/build

# executable must be present in the current directory! #
export rundir=`pwd`
echo "rundir=$rundir"

#lcg-cp lfn://grid/voce/ilias/dirac.x file://$rundir/dirac.x  || { echo "Exiting...no dirac.x found..." ; exit 2; }
lcg-cp lfn://grid/compchem/ilias/dirac.x file://$rundir/dirac.x  || { echo "Exiting...no dirac.x found..." ; exit 2; }

##########################################
##  make it executable and finally launch
##########################################
chmod 755 dirac.x
export DIRWRK=84000000
./dirac.x
