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
echo "Zero field..."
time $PAM_PATH/pam -rmscr -rmwrk -mw 10 $MPIARG -utcmo \
   HF.p-ccpvdz.mol zerofield.inp
echo "No relaxation..."
time $PAM_PATH/pam -rmscr -rmwrk -mw 10 $MPIARG -incmo \
   HF.p-ccpvdz.mol ccsd.1
time $PAM_PATH/pam -rmscr -rmwrk -mw 10 $MPIARG -incmo \
   HF.p-ccpvdz.mol ccsd.m1
echo "With relatixation"
time $PAM_PATH/pam -rmscr -rmwrk -mw 10 $MPIARG -incmo -utcmo \
   HF.p-ccpvdz.mol scf.ccsd.1
time $PAM_PATH/pam -rmscr -rmwrk -mw 10 $MPIARG -incmo \
   HF.p-ccpvdz.mol scf.ccsd.m1
rm DFCOEF
