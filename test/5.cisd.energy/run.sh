#!/bin/sh
#
# Set up calculation
#

MOL=F;
export MOL;
FILTER=filter
export FILTER;

###############################
echo "DC cisd run"
MENU=cisd;
export MENU;
../run.sh
###############################
echo "BSS cisd run"
MENU=F.bss.cisd;
export MENU;
../run.sh
###############################
echo "BSS+MFSSO cisd run"
MENU=F.bss+mfsso.cisd;
export MENU;
../run.sh
