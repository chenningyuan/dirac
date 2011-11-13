#!/bin/sh
#
# Set up calculation
#  
MOL=Hg;
export MOL;
#
MENU=std;
export MENU;
../run.sh
#
MENU=free;
export MENU;
../run.sh
