#!/bin/sh
#
# Set up calculation
#
# MI: There are too many two-component Hamiltonians,
# resort tests into 3 main groups to reduce the computer time
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
echo "BSS+MFSSO cisd run"
MENU=F.bss+mfsso.cisd;
export MENU;
../run.sh
###############################
echo "BSSsfESO1+MFSSO cisd run"
MENU=F.bss_sfEso1+mfsso.cisd;
export MENU;
../run.sh
###############################
if [ "$RUNTYPE" = "short" ]; then
   echo
   echo "Completed short test"
   echo "<-----"
   exit
fi
###############################
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
echo "DC2BSS_RKB(DF) cisd run"
MENU=F.dc2bss_rkb_df.cisd;
export MENU;
../run.sh
###############################
echo "DKH2sfESO1+MFSSO cisd run"
MENU=F.dkh2_sfEso1+mfsso.cisd;
export MENU;
../run.sh
###############################

###############################
if [ "$RUNTYPE" = "medium" ]; then
   echo
   echo "Completed medium test"
   echo "<-----"
   exit
fi
###############################
echo "BSS_RKB+MFSSO cisd run"
MENU=F.bss_rkb+mfsso.cisd;
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
if [ "$RUNTYPE" = "long" ]; then
   echo
   echo "Completed long test"
   echo "<-----"
   exit
fi
###############################

