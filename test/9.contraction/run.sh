#!/bin/sh
#
# Set up calculation
#
MENU=dhf;
MOL=H2O.large+small;
export MENU;
export MOL;
../run.sh
MENU=dhf;
MOL=H2O.usingfiles;
export MENU;
export MOL;
../run.sh -copy O.ccpvdz.bas -copy H.ccpvdz.bas
