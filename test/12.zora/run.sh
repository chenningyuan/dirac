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
time $PAM_PATH/pam -rmscr -rmwrk $MPIARG -mw 10 U.dyall.2p 1el.dirac
time $PAM_PATH/pam -rmscr -rmwrk $MPIARG -mw 10 U.dyall.2p 1el.diracsf
time $PAM_PATH/pam -rmscr -rmwrk $MPIARG -mw 10 U.dyall.2p 1el.zora
time $PAM_PATH/pam -rmscr -rmwrk $MPIARG -mw 10 U.dyall.2p 1el.zorasf
time $PAM_PATH/pam -rmscr -rmwrk $MPIARG -mw 10 U.dyall.2p 1el.zorasc
time $PAM_PATH/pam -rmscr -rmwrk $MPIARG -mw 10 U.dyall.2p 1el.zorascsf
time $PAM_PATH/pam -rmscr -rmwrk $MPIARG -mw 10 Xe.Faas nonrel
time $PAM_PATH/pam -rmscr -rmwrk $MPIARG -mw 10 Xe.Faas zora
time $PAM_PATH/pam -rmscr -rmwrk $MPIARG -mw 10 Xe.Faas zora4
time $PAM_PATH/pam -rmscr -rmwrk $MPIARG -mw 10 Xe.Faas dirac
time $PAM_PATH/pam -rmscr -rmwrk $MPIARG -mw 10 Xe.Faas zorasf
time $PAM_PATH/pam -rmscr -rmwrk $MPIARG -mw 10 Xe.Faas zora4sf
time $PAM_PATH/pam -rmscr -rmwrk $MPIARG -mw 10 Xe.Faas diracsf
time $PAM_PATH/pam -rmscr -rmwrk $MPIARG -mw 10 Xe.Faas zorasc
time $PAM_PATH/pam -rmscr -rmwrk $MPIARG -mw 10 Xe.Faas zora4sc
time $PAM_PATH/pam -rmscr -rmwrk $MPIARG -mw 10 Xe.Faas zorascsf
time $PAM_PATH/pam -rmscr -rmwrk $MPIARG -mw 10 Xe.Faas zora4scsf
