#!/bin/sh
#
# Set up calculation
#  
MOL=H2O2; export MOL;
MENU=shielding.dc; export MENU;
../run.sh
MENU=shielding.ll; export MENU;
../run.sh
MENU=spinspin.dc; export MENU;
../run.sh
MENU=spinspin.ll; export MENU;
../run.sh
