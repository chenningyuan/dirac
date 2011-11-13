#!/bin/sh
#
# Set up calculation
#
MOL=H2O;
export MOL;
#######################
echo "Dirac-Coulomb Hamiltonian"
MENU=rpa;
export MENU;
../run.sh
#######################
#echo "BSS+MFSSO  Hamiltonian"
#MENU=H2O.bss+mfsso.rpa;
#export MENU;
#../run.sh
#######################
#echo "BSS_rkb+MFSSO  Hamiltonian"
#MENU=H2O.bss_rkb+mfsso.rpa;
#export MENU;
#../run.sh
#######################
echo "BSS  Hamiltonian"
MENU=H2O.bss.rpa;
export MENU;
../run.sh
#######################
