#!/bin/sh
#
# Set up calculation
#
MOL=Z80H.lsym.dir.mol;
export MOL;
###############################
echo "bss.ccsd run"
FILTER=filter
MENU=HF.bss.cc;
export FILTER;
export MENU;
../run.sh
###############################
echo "bss_sfb.ccsd-ci run"
FILTER=filter_cc+ci
MENU=HF.bss_sfb.cc_cisd;
export FILTER;
export MENU;
../run.sh
###############################
echo "dk2_sfb.ccsd-ci run"
MENU=HF.dk2_sfb.pnuc.cc_cisd;
export MENU;
../run.sh
