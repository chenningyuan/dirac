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
echo "BSSsfBSO1 Hamiltonian"
MENU=bss_sfBso1.ci;
export MENU;
../run.sh
######################################"
echo "BSSsfESO1 Hamiltonian"
MENU=bss_sfEso1.ci;
export MENU;
../run.sh
######################################"
echo "BSS_rkb Hamiltonian"
MENU=bss_rkb.ci;
export MENU;
../run.sh
######################################"
echo "BSS+MFSSO Hamiltonian"
MENU=bss+mfsso.ci;
export MENU;
../run.sh
######################################"
echo "BSSsfBSO1+MFSSO Hamiltonian"
MENU=bss_sfBso1+mfsso.ci;
export MENU;
../run.sh
######################################"
echo "BSSsfESO1+MFSSO Hamiltonian"
MENU=bss_sfEso1+mfsso.ci;
export MENU;
../run.sh
######################################"
echo "BSS_rkb+MFSSO Hamiltonian"
MENU=bss_rkb+mfsso.ci;
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
echo "DKH2sfBSO1+MFSSO Hamiltonian"
MENU=dkh2sfBso1+mfsso.ci;
export MENU;
../run.sh
######################################"
echo "DKH2sfESO1+MFSSO Hamiltonian"
MENU=dkh2sfEso1+mfsso.ci;
export MENU;
../run.sh
######################################"
echo "DKH2sfBSO1 Hamiltonian"
MENU=dkh2sfBso1.ci;
export MENU;
../run.sh
######################################"
echo "DKH2sfESO1 Hamiltonian"
MENU=dkh2sfEso1.ci;
export MENU;
../run.sh
######################################"

