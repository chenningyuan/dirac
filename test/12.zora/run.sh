#!/bin/sh
#
# Set up calculation
#
MOL=U.dyall.2p;
export MOL;
MENU=1el.dirac;
export MENU;
../run.sh
MENU=1el.diracsf;
export MENU;
../run.sh
MENU=1el.zora;
export MENU;
../run.sh
MENU=1el.zorasf;
export MENU;
../run.sh
MENU=1el.zorasc;
export MENU;
../run.sh
MENU=1el.zorascsf;
export MENU;
../run.sh
