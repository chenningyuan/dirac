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
../run.sh
#############
echo "DC Hamiltonian, IH FSCCSD"
MENU=fsccsd_IH
export MENU;
../run.sh
#############
#REMOVE_BEGIN
#echo "BSS Hamiltonian, IH FSCCSD"
#MENU=bss.fsccsd_IH;
#export MENU;
#../run.sh
#REMOVE_END
#############
#REMOVE_BEGIN
echo "BSS+MFSSO Hamiltonian, IH FSCCSD"
MENU=bss+mfsso.fsccsd_IH;
export MENU;
../run.sh
#REMOVE_END
#############

