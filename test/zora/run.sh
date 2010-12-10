#!/bin/sh
#
# Set up calculation
#
MOL=U.dyall.2p.mol;
export MOL;
MENU=1el.dirac;
export MENU;
../SCRIPTS/run.sh
MENU=1el.diracsf;
export MENU;
../SCRIPTS/run.sh
MENU=1el.zora;
export MENU;
../SCRIPTS/run.sh
MENU=1el.zorasf;
export MENU;
../SCRIPTS/run.sh
MENU=1el.zorasc;
export MENU;
../SCRIPTS/run.sh
MENU=1el.zorascsf;
export MENU;
../SCRIPTS/run.sh
