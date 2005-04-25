#!/bin/sh
#
# Set up calculation
#
MOL=F;
export MOL;
######################################"
echo "DC Hamiltonian"
MENU=ci;
export MENU;
../run.sh
######################################"
echo "BSS Hamiltonian"
MENU=bss.ci;
export MENU;
../run.sh
######################################"
echo "BSS+MFSSO Hamiltonian"
MENU=bss+mfsso.ci;
export MENU;
../run.sh
######################################"
echo "DKH2 Hamiltonian"
MENU=dkh2.ci;
export MENU;
../run.sh
######################################"
echo "DKH2+MFSSO Hamiltonian"
MENU=dkh2+mfsso.ci;
export MENU;
../run.sh
######################################"
echo "DKH2+SO1 Hamiltonian"
MENU=dkh2so1.ci;
export MENU;
../run.sh
######################################"

