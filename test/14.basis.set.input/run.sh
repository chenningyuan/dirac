#!/bin/sh
#
# Set up calculation
#  
MENU=ch4; export MENU;
MOL=ch4_temp; export MOL;
../run.sh
MOL=ch4_explicit; export MOL;
../run.sh
MOL=ch4_family; export MOL;
../run.sh
MOL=ch4_dualfamily; export MOL;
../run.sh
MOL=ch4_basis; export MOL;
../run.sh -basdir `pwd`
MOL=ch4_molfbas; export MOL;
../run.sh -basdir `pwd`:`pwd`/../../basis/MOLFDIR/accpvtz
