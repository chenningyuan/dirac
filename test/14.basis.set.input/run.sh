#!/bin/sh
#
# Set up calculation
#  
MENU=ch4 export MENU;
MOL=ch4_temp; export MOL;
../run.sh
MOL=ch4_explicit; export MOL;
../run.sh
MOL=ch4_family; export MOL;
../run.sh
MOL=dualfamily; export MOL;
../run.sh
#time $PAM_PATH/pam -rmscr -rmwrk -mw 3 $MPIARG -basdir `pwd` ch4_basis ch4
#time $PAM_PATH/pam -rmscr -rmwrk -mw 3 $MPIARG -basdir `pwd` \
#   -basdir `pwd`/../../basis/MOLFDIR/accpvtz ch4_molfbas ch4
