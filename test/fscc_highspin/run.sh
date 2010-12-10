#!/bin/sh
#
# Set up calculation
#  
MOL=O2;
export MOL;
#############
echo "COSCI"
MENU=hf;
export MENU;
../SCRIPTS/run.sh
#############
echo "Unrestricted CCSD with NELEC"
MENU=uccsd;
export MENU;
../SCRIPTS/run.sh
#############
echo "Unrestricted CCSD with NELEC_F1, NELEC_F2"
MENU=uccsd_linsym;
export MENU;
../SCRIPTS/run.sh
#############
echo "Fock space CCSD"
MENU=fsccsd;
export MENU;
../SCRIPTS/run.sh
