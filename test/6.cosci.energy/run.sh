#!/bin/ksh
#
# Set up calculation
#
MENU=ci
MOL=F
if [ -z "$PAM_PATH" ]; then
   PAM_PATH=`which pam`
   if [ -n "`echo $PAM_PATH | grep 'Command not found'`" ]; then
      echo "pam not found. Please either set the enviroment variable"
      echo "PAM_PATH or put pam in your path."
      exit 1
   fi
   PAM_PATH=${PAM_PATH%/*}
fi
if [ -z "$DIRAC" ]; then
   DIRAC=`which dirac.x`
fi
if [ -z "$TAG" ]; then
   TAG=out
fi
NEWOUT=${MENU}_${MOL}.${TAG}
OLDOUT=${MENU}_${MOL}.out 
#
# Run calculation
#
time $PAM_PATH/pam -dirac $DIRAC -run $TAG -rmwrk -rmscr $MPIARG $MOL $MENU
#
# Generate reports
#
../crep  result/$OLDOUT result/filter >result/$OLDOUT.rep
../crep  $NEWOUT result/filter >$NEWOUT.rep
#
# Generate differences
#
../bdrep $NEWOUT.rep result/$OLDOUT.rep >${NEWOUT}_dif.txt
../drep  $NEWOUT.rep result/$OLDOUT.rep >${NEWOUT}_dif.htm
#
# Add on to grand total of differences
#
../hdrep  $NEWOUT.rep result/$OLDOUT.rep >>../${TAG}_dif.htm
cat $NEWOUT_dif.txt >>../${TAG}_dif.txt
