#!/bin/sh
#
# Set up calculation
#  
MOL=Hg;
export MOL;
#
MENU=std;
export MENU;
../SCRIPTS/run.sh
#
MENU=free;
export MENU;
../SCRIPTS/run.sh
