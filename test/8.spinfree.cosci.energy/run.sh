#!/bin/sh
#
# Set up calculation
#

#######################################
echo "DCspinfree calculations"
MENU=ci;
MOL=F;
export MENU;
export MOL;
../run.sh
#######################################
echo "BSS_rkb_sfB calculations"
MENU=bss_rkb_sfB.ci;
MOL=F;
export MENU;
export MOL;
../run.sh

###############################
if [ "$RUNTYPE" = "short" ]; then
   echo
   echo "Completed short test"
   echo "<-----"
   exit
fi
###############################

#######################################
echo "BSS_rkb_sfE calculations"
MENU=bss_rkb_sfE.ci;
MOL=F;
export MENU;
export MOL;
../run.sh
#######################################
echo "BSS_sfB calculations"
MENU=bss_sfB.ci;
MOL=F;
export MENU;
export MOL;
../run.sh
###############################
if [ "$RUNTYPE" = "medium" ]; then
   echo
   echo "Completed medium test"
   echo "<-----"
   exit
fi
#skipping outdated BSS tests : update output if you want them...
exit
###############################
#######################################
echo "BSS_sfE calculations"
MENU=bss_sfE.ci;
MOL=F;
export MENU;
export MOL;
../run.sh
#######################################
echo "DKH2_sfB calculations"
MENU=dkh2sfB.ci;
MOL=F.lsym.mol;
export MENU;
export MOL;
../run.sh
#######################################

