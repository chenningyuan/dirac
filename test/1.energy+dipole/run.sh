#!/bin/ksh
if [ -z "$PAM_PATH" ]; then
   PAM_PATH=`/usr/bsd/which pam`
   if [ -n "`echo $PAM_PATH | grep 'Command not found'`" ]; then
      echo "pam not found. Please either set the enviroment variable"
      echo "PAM_PATH or put pam in your path."
      exit 1
   fi
   PAM_PATH=${PAM_PATH%/*}
fi
time $PAM_PATH/pam -rmscr -rmwrk -mw 3 $MPIARG H2O ccsd
