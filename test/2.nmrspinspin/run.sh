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
../run.sh
#################################
echo "BSS Hamiltonian"
MENU=bss.spsp;
export MENU;
../run.sh
#################################
#################################
