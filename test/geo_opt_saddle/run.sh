#!/bin/sh
#
# Set up calculation
#  
MOL=ethane_1; export MOL;
MENU=minimum; export MENU;
../SCRIPTS/run.sh
MENU=saddle; export MENU;
../SCRIPTS/run.sh
MOL=ethane_2; export MOL;
MENU=findre; export MENU;
../SCRIPTS/run.sh
MENU=constrained; export MENU;
../SCRIPTS/run.sh

