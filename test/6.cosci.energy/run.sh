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

MOL=F;
export MOL;
######################################"
echo "DC Hamiltonian"
MENU=ci;
export MENU;
../run.sh
######################################"
echo "BSS Hamiltonian"
MENU=bss.ci;
export MENU;
../run.sh
######################################"
###############################
if [ "$RUNTYPE" = "short" ]; then
   echo
   echo "Completed short test"
   echo "<-----"
   exit
fi
###############################
######################################"
echo "BSSsfBSO1 Hamiltonian"
MENU=bss_sfBso1.ci;
export MENU;
../run.sh
######################################"
echo "BSSsfESO1 Hamiltonian"
MENU=bss_sfEso1.ci;
export MENU;
../run.sh
######################################"
echo "BSS_rkb Hamiltonian"
MENU=bss_rkb.ci;
export MENU;
../run.sh
######################################"
echo "BSS+MFSSO Hamiltonian"
MENU=bss+mfsso.ci;
export MENU;
../run.sh
###############################
if [ "$RUNTYPE" = "medium" ]; then
   echo
   echo "Completed medium test"
   echo "<-----"
   exit
fi
###############################
######################################"
echo "BSSsfBSO1+MFSSO Hamiltonian"
MENU=bss_sfBso1+mfsso.ci;
export MENU;
../run.sh
######################################"
echo "BSSsfESO1+MFSSO Hamiltonian"
MENU=bss_sfEso1+mfsso.ci;
export MENU;
../run.sh
######################################"
echo "BSS_rkb+MFSSO Hamiltonian"
MENU=bss_rkb+mfsso.ci;
export MENU;
../run.sh
######################################"
echo "DKH2 Hamiltonian"
MENU=dkh2.ci;
export MENU;
../run.sh
######################################"
echo "DKH2+MFSSO Hamiltonian"
MENU=dkh2+mfsso.ci;
export MENU;
../run.sh
######################################"
echo "DKH2sfBSO1+MFSSO Hamiltonian"
MENU=dkh2sfBso1+mfsso.ci;
export MENU;
../run.sh
######################################"
echo "DKH2sfESO1+MFSSO Hamiltonian"
MENU=dkh2sfEso1+mfsso.ci;
export MENU;
../run.sh
######################################"
echo "DKH2sfBSO1 Hamiltonian"
MENU=dkh2sfBso1.ci;
export MENU;
../run.sh
######################################"
echo "DKH2sfESO1 Hamiltonian"
MENU=dkh2sfEso1.ci;
export MENU;
../run.sh
######################################"
### Heavy F analogue
######################################"
MOL=Z86_Flike;
export MOL;
######################################"
echo "DC Hamiltonian, Rn(77+)"
MENU=ci;
export MENU;
../run.sh
######################################"
echo "BSS Hamiltonian, Rn(77+)"
MENU=bss.ci;
export MENU;
../run.sh
######################################"
echo "DKH2 Hamiltonian, Rn(77+)"
MENU=dkh2.ci;
export MENU;
../run.sh
######################################"
echo "BSS_RKB (TROND) Hamiltonian, Rn(77+)"
MENU=bss_rkb.ci;
export MENU;
../run.sh
######################################"
echo "DKH2+MFSSO Hamiltonian, Rn(77+)"
MENU=dkh2+mfsso_Flike.ci;
export MENU;
../run.sh
######################################"
echo "DKH2sfBSO1 Hamiltonian, Rn(77+)"
MENU=dkh2sfBso1.ci;
export MENU;
../run.sh
######################################"
echo "DKH2sfBSO1+MFSSO Hamiltonian, Rn(77+)"
MENU=dkh2sfBso1+mfsso.ci_Flike;
export MENU;
../run.sh
######################################"
echo "BSS+MFSSO Hamiltonian, Rn(77+)"
MENU=bss+mfsso.ci_Flike;
export MENU;
../run.sh
######################################"
###############################
if [ "$RUNTYPE" = "long" ]; then
   echo
   echo "Completed long test"
   echo "<-----"
   exit
fi
###############################
