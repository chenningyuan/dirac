#!/bin/ksh
#
# Generic place for the log-file
#
LOG=../logs
#
# Set the default arguments for the pam script
#
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
if [ "$MOL" = "$MENU" ]; then
   NEWOUT=${MENU}.${TAG}
   OLDOUT=${MENU}.out 
else
   NEWOUT=${MENU}_${MOL}.${TAG}
   OLDOUT=${MENU}_${MOL}.out 
fi
#
# Run calculation (pass on additional arguments for pam that are supplied by the caller)
#
EXIT=0 
$PAM_PATH/pam -dirac $DIRAC -run $TAG -rmwrk -rmscr $MPIARG $* $MOL $MENU || EXIT=1
if [ $EXIT -eq 0 ]; then
#
# Succesful run : check for deviations with the standard output.
#
# Extract report from the output using the filter file
#
if [ -z "$FILTER" ]; then
   FILTER=filter
fi
../crep  result/$OLDOUT result/$FILTER >result/$OLDOUT.rep
../crep  $NEWOUT result/$FILTER >$NEWOUT.rep
#
# Generate differences
#
../bdrep $NEWOUT.rep result/$OLDOUT.rep >${NEWOUT}_dif.txt
../drep  $NEWOUT.rep result/$OLDOUT.rep >${NEWOUT}_dif.htm
#
# Add on to grand total of differences
#
../hdrep  $NEWOUT.rep result/$OLDOUT.rep >>${LOG}/${TAG}_dif.htm
cat ${NEWOUT}_dif.txt                    >>${LOG}/${TAG}_dif.txt
#
else
#
# Run crashed : write report.
#
echo " <HR>"                                                    >>${LOG}/${TAG}_dif.htm
echo "<FONT COLOR="#FF0000"> `pwd` : $MOL $MENU failed </FONT>" >>${LOG}/${TAG}_dif.htm
echo " <HR>"                                                    >>${LOG}/${TAG}_dif.htm
echo " `pwd` : $MOL $MENU failed"                               >>${LOG}/${TAG}_dif.txt
fi
#
