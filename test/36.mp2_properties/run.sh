#!/bin/sh
#
# Set up calculation
#
MOL=H2O2;
export MOL;
###############################
echo "MP2 first order properties..."
MENU=mp2_gradient;
FILTER=filter;
export MENU;
export FILTER;
../run.sh
###############################
