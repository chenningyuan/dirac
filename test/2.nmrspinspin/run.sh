#!/bin/sh
#
# Set up calculation
#
#################################
MOL=H2O;
export MOL;
#################################
echo "DC Hamiltonian"
MENU=spsp;
export MENU;
../run.sh -mw 100
#################################
echo "BSS Hamiltonian"
MENU=bss.spsp;
export MENU;
../run.sh -mw 100
#################################
#REMOVE_BEGIN
echo "BSS+MFSSO Hamiltonian"
MENU=bss+mfsso.spsp;
export MENU;
../run.sh -mw 100
#################################
echo "BSS_rkb+MFSSO Hamiltonian"
MENU=bss_rkb+mfsso.spsp;
export MENU;
../run.sh -mw 100
#REMOVE_END
#################################
