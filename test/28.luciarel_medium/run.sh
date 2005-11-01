#!/bin/sh
#
# Set up calculation
#
MOL=F;
export MOL;
#  DC Hamiltonian ...
echo "DC ci run"
MENU=F;
export MENU;
../run.sh
### BSS+MFSSO Hamiltonian
echo "BSS+MFSSO ci run"
MENU=F.bss+mfsso;
export MENU;
../run.sh

