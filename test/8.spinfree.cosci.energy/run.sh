#!/bin/sh
#
# Set up calculation
#

if [ ! -z "$RUNTYPE" ]; then
   case "${RUNTYPE}" in
      "short"|"medium"|"long")
#        MENU=cc
         ;;
      *)
      echo "we do not run the benchmark with runtype $RUNTYPE" 
         exit
         ;;
   esac
else
   echo "we do not run this test if the runtype is not specified !" 
   exit
fi

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
MOL=F.lsym;
export MENU;
export MOL;
../run.sh
#######################################

