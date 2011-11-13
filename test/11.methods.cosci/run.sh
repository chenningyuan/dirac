#!/bin/sh
#
# Set up calculation
#
MOL=O;
export MOL;
###############################
MENU=ci.1;
export MENU;
../run.sh
###############################
MENU=ci.2;
export MENU;
../run.sh -utcmo
###############################
MENU=ci.3;
export MENU;
../run.sh -incmo
###############################
MENU=ci.4;
export MENU;
../run.sh
###############################
rm DFCOEF
