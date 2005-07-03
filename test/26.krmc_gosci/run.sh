#!/bin/sh
#
# Set up calculation
#
MOL=Be;
export MOL;
##### DC ######
echo "Dirac-Coulomb KR-MC-SCF calculations"
MENU=Be;
export MENU;
DIRWRK=30000000; export DIRWRK
../run.sh
#######################################
echo "BSS+MFSSO KR-MC-SCF calculations"
MENU=Be.bss+mfsso.krmc;
export MENU;
DIRWRK=30000000; export DIRWRK
../run.sh
#######################################

