#!/bin/sh
#
# Set up calculation
#  
MENU=std;
MOL=Hg;
export MENU;
export MOL;
../run.sh
MENU=free;
export MENU;
../run.sh
