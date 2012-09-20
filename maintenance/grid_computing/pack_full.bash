#!/bin/bash
#
# simple script for preparing Dirac tar-ball for grid computing
#


# position of your DIRAC trunk
DIRAC=/home/ilias/qch_work/software/dirac_git/trunk_current

### update & pack ###

cd $DIRAC

#echo -e "\n --- obtaining the fresh trunk ---"
#git pull
#echo -e "\n --- last git revision ---"
#git log -n 1

#------------------------------------------------------------------------------
# make one tgz-file with all Dirac stuff for testing & running on a grid CE
#------------------------------------------------------------------------------

# problems with "--exclude test/*/*.out"  - so remove them
#/bin/rm test/*/*.out*
#/bin/rm -r test/*/compare

if [[ -e $DIRAC/maintenance/dirac_current.tgz ]];  then
   echo -e "\n removing old dirac_current.tgz file..."
   /bin/rm dirac_current.tgz
fi


# pack all necessary stuff

echo -e "\n now packing... "

tar cvzf maintenance/grid_testing/dirac_current.tgz  --exclude="maintenance/grid_testing/dirac_current.tgz test/*/compare  test/*/*.out"  *

echo -e "\n --- package dirac_current.tgz is ready, you can copy it onto chosen grid SE ---"
echo    "-----------------------------------------------------------------------------------"

echo -e "\n DIRAC package is of git revision :"
git log -n 1

exit 0
