#!/bin/sh
#
# Set up calculation
#
MENU=levy-leblond;
MOL=H2O;
export MENU;
export MOL;
../run.sh
MENU=nmr_LinK_Jengine
MOL=H2O_nosym
export MENU;
export MOL;
../run.sh
MENU=nmr_hermit
MOL=H2O_nosym
export MENU;
export MOL;
../run.sh
MENU=nmr_LinK_df-Jengine
MOL=H2O_df
export MENU;
export MOL;
../run.sh
