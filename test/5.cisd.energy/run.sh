#!/bin/sh
#
# Set up calculation
#

MOL=F;
export MOL;
FILTER=filter
export FILTER;

###############################
echo "DC cisd run"
MENU=cisd;
export MENU;
../run.sh
###############################
echo "BSS cisd run"
MENU=F.bss.cisd;
export MENU;
../run.sh
###############################
echo "BSS_RKB cisd run"
MENU=F.bss_rkb.cisd;
export MENU;
../run.sh
###############################
echo "BSS+MFSSO cisd run"
MENU=F.bss+mfsso.cisd;
export MENU;
../run.sh
###############################
echo "BSS_RKB+MFSSO cisd run"
MENU=F.bss_rkb+mfsso.cisd;
export MENU;
../run.sh
###############################
echo "BSSsfESO1+MFSSO cisd run"
MENU=F.bss_sfEso1+mfsso.cisd;
export MENU;
../run.sh
###############################
echo "BSSsfESO1 cisd run"
MENU=F.bss_sfEso1.cisd;
export MENU;
../run.sh
###############################
echo "BSSsfBSO1 cisd run"
MENU=F.bss_sfBso1.cisd;
export MENU;
../run.sh
###############################
echo "DKH2sfBSO1 cisd run"
MENU=F.dkh2_sfBso1.cisd;
export MENU;
../run.sh
###############################
echo "DKH2sfESO1 cisd run"
MENU=F.dkh2_sfEso1.cisd;
export MENU;
../run.sh
###############################
echo "DKH2sfESO1+MFSSO cisd run"
MENU=F.dkh2_sfEso1+mfsso.cisd;
export MENU;
../run.sh
###############################
