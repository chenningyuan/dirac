#!/bin/sh
#
# Set up calculation
#
MOL=H2O2;
export MOL;
FILTER=filter;
export FILTER;
###############################
echo "MP2 first order properties...DC Hamiltonian"
MENU=mp2_gradient;
export MENU;
../run.sh
###############################
