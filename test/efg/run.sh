#!/bin/sh
#
# Set up calculation
#
MOL=HF.p-ccpvdz.mol;
export MOL;
###############################
echo "Zero field...DC"
MENU=zerofield;
export MENU;
../SCRIPTS/run.sh --mw=140  --aw=160 --outcmo
###############################
echo "Field strength 0.001 ... DC"
MENU=scf.ccsd.1;
export MENU;
../SCRIPTS/run.sh --mw=140 --aw=160 --incmo
###############################
echo "Field strength -0.001 ... DC"
MENU=scf.ccsd.m1;
export MENU;
../SCRIPTS/run.sh --mw=140  --aw=160 --incmo
###############################
rm DFCOEF
###############################
###############################
#<REMOVE_BEGIN>
echo "Zero field...BSS+MFSSO "
MENU=bss+mfsso.zerofield;
export MENU;
../SCRIPTS/run.sh --outcmo
#<REMOVE_END>
###############################
echo "Zero field...BSS "
MENU=bss.zerofield;
export MENU;
../SCRIPTS/run.sh --incmo --outcmo
###############################
echo "Field strength 0.001 ... BSS"
MENU=bss.scf.ccsd.1;
export MENU;
../SCRIPTS/run.sh --incmo --mw=140 --aw=160
###############################
echo "Field strength -0.001 ... BSS"
MENU=bss.scf.ccsd.m1;
export MENU;
../SCRIPTS/run.sh --incmo --mw=140 --aw=160
###############################
rm DFCOEF

