#!/bin/sh
#
# Set up calculation
#
MENU=dhf;
MOL=H2O.large+small.mol;
export MENU;
export MOL;
../SCRIPTS/run.sh
MENU=dhf;
MOL=H2O.usingfiles.mol;
export MENU;
export MOL;
../SCRIPTS/run.sh --copy="O.ccpvdz.bas H.ccpvdz.bas"
