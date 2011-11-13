#!/bin/sh
#
# Set up calculation
#  
MOL=ethane_1; export MOL;
MENU=minimum; export MENU;
../run.sh
MENU=saddle; export MENU;
../run.sh
MOL=ethane_2; export MOL;
MENU=findre; export MENU;
../run.sh
MENU=constrained; export MENU;
../run.sh

