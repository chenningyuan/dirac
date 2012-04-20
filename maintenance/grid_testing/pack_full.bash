#!/bin/bash
#
# simple script for preparing Dirac tar-ball for grid computing
#


DIRAC=/home/ilias/qch_work/software/dirac_git/trunk_current

### update & pack ###

cd $DIRAC

#echo -e "\n --- obtaining the fresh trunk ---"
#git pull

#echo -e "\n --- last git revision ---"
#git log -n 1

#--------------------------------------------------------------
# make one tgz-file with all Dirac stuff for testing & running 
#--------------------------------------------------------------

# problems with "--exclude test/*/*.out"  - so remove them
/bin/rm test/*/*.out*
/bin/rm test/*/compare

if [[ -e $DIRAC/maintenance/dirac_current.tgz ]];  then
   echo -e "\n removing old dirac_current.tgz..."
   /bin/rm dirac_current.tgz
fi

echo -e "\n now packing... "

# pack all necessary stuff
tar cvzf maintenance/grid_testing/dirac_current.tgz  ----exclude maintenance/grid_testing/dirac_current.tgz  *
#mv   $DIRAC/dirac_current.tgz   maintenance/grid_testing/.

echo -e "\n --- package dirac_current.tgz ready, you can copy it onto chosen grid SE ---"

echo -e "\n --- and the last git revision ---"
git  log -n 1

exit 0
