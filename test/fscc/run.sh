#!/bin/sh
#
# Set up calculation
#  
MOL=Mg;
export MOL;
#############
echo "DC Hamiltonian, FSCCSD"
MENU=fsccsd
export MENU;
../run.sh -mw 100
#############
echo "BSS+MFSSO Hamiltonian, IH FSCCSD"
MENU=bss+mfsso.fsccsd_IH
export MENU;
../run.sh -mw 110
#############
echo "DC Hamiltonian, IH FSCCSD"
MENU=fsccsd_IH
export MENU;
../run.sh -mw 110
#############
