#!/bin/ksh
if [ -z "$PAM_PATH" ]; then
   PAM_PATH=`which pam`
   if [ -n "`echo $PAM_PATH | grep 'Command not found'`" ]; then
      echo "pam not found. Please either set the enviroment variable"
      echo "PAM_PATH or put pam in your path."
      exit 1
   fi
   PAM_PATH=${PAM_PATH%/*}
fi
time $PAM_PATH/pam -rmscr -rmwrk -basdir `pwd` -mw 4 $MPIARG ethane_1 minimum
time $PAM_PATH/pam -rmscr -rmwrk -basdir `pwd` -mw 4 $MPIARG ethane_1 saddle
time $PAM_PATH/pam -rmscr -rmwrk -basdir `pwd` -mw 4 $MPIARG ethane_2 findre
time $PAM_PATH/pam -rmscr -rmwrk -basdir `pwd` -mw 4 $MPIARG ethane_2 constrained
exit 0

