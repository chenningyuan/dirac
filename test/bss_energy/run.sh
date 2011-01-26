#!/bin/sh
#
# Set up calculation
#
MOL=Z80H.lsym.dir.mol;
export MOL;
###############################
echo "bss.ccsd run"
MENU=HF.bss.cc;
export MENU;
../SCRIPTS/run.sh --mw=80
###############################
echo "bss_sfb.ccsd-ci run"
MENU=HF.bss_sfb.cc_cisd;
export MENU;
../SCRIPTS/run.sh
###############################
echo "dk2_sfb.ccsd-ci run"
MENU=HF.dk2_sfb.pnuc.cc_cisd;
export MENU;
../SCRIPTS/run.sh
