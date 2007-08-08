#!/bin/sh
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
#
# Run calculation (pass on additional arguments for pam that are supplied by the caller)
#
EXIT=0 
$PAM_PATH/pam -dirac $DIRAC $NICE -run $TAG -rmwrk -clean $MPIARG $* $MOL $MENU || EXIT=1
#
# Finished : either regular or with a crash. Make two branches.
#
if [ $EXIT -eq 0 ]; then
#
#  Succesful run : check for deviations with the standard output.
#
MENU=${MENU%\.inp}
MOL=${MOL%\.mol}
MOL=${MOL%\.xyz}
if [ "$MOL" = "$MENU" ]; then
   NEWOUT=${MENU}.${TAG}
   OLDOUT=${MENU}.out 
else
   NEWOUT=${MENU}_${MOL}.${TAG}
   OLDOUT=${MENU}_${MOL}.out 
fi
#
#  Extract report from the output using the filter file
#
  if [ -z "$FILTER" ]; then
     FILTER=filter
  fi
  ../crep  result/$OLDOUT result/$FILTER >result/$OLDOUT.rep
  ../crep  $NEWOUT result/$FILTER >$NEWOUT.rep
#
#  Generate file with the observed discrepancies with the standard output that is placed
#  in the directory of the test. This is done both in html format and in txt format.
#
  ../bdrep $NEWOUT.rep result/$OLDOUT.rep >${NEWOUT}_dif.txt
  ../drep  $NEWOUT.rep result/$OLDOUT.rep >${NEWOUT}_dif.htm
#
#  Now we take this detailed output and add it on to the grand total of differences.
#  The html output is added to a file that contains ALL the selected output as defined by
#  the filter files. Observed discrepancies ( possible errors) are colour-coded.
#
  ../hdrep  $NEWOUT.rep result/$OLDOUT.rep >>${LOG}/${TAG}_dif.htm
#
#  We also write a SINGLE LINE of output ONLY IF a test fails to total txt file. 
#  The user should then check the full output in the directory of the specific run. 
#  This compact output makes the visual inspection of the output very easy and allows for 
#  regular daily use (e.g. via a crontab script).
#
   if [ -s ${NEWOUT}_dif.txt ]; then
      echo " Test `pwd` $MOL $MENU failed : output not consistent with reference output" >>${LOG}/${TAG}_dif.txt
   else
#
#     This strange construction (cat an empty file) is needed to make sure that the log file
#     exists to enable similar checking (using if -e and then if -s) in the testlast script.
#
      cat ${NEWOUT}_dif.txt >>${LOG}/${TAG}_dif.txt
   fi
#
else
#
# Run crashed : write report.
#
# Write only write a report in the logs directory.
#
  echo " <HR>"                                                    >>${LOG}/${TAG}_dif.htm
  echo "<FONT COLOR="#FF0000"> `pwd` : $MOL $MENU failed </FONT>" >>${LOG}/${TAG}_dif.htm
  echo " <HR>"                                                    >>${LOG}/${TAG}_dif.htm
  echo " Test `pwd` $MOL $MENU failed : runtime error (DIRAC crashed)" >>${LOG}/${TAG}_dif.txt
#
fi
#
