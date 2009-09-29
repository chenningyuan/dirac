#!/bin/sh
#
# Set up calculation
#
MOL=HF.p-ccpvdz.mol;
export MOL;
###############################
echo "Zero field...DC"
MENU=zerofield;
FILTER=filter_zerofield;
export MENU;
export FILTER;
../run.sh -mw 140 -utcmo
###############################
echo "Field strength 0.001 ... DC"
MENU=scf.ccsd.1;
export MENU;
../run.sh -mw 140 -incmo
###############################
echo "Field strength -0.001 ... DC"
MENU=scf.ccsd.m1;
export MENU;
../run.sh -mw 140 -incmo
###############################
rm DFCOEF
###############################
###############################
#<REMOVE_BEGIN>
echo "Zero field...BSS+MFSSO "
MENU=bss+mfsso.zerofield;
FILTER=filter_zerofield;
export MENU;
export FILTER;
../run.sh -utcmo
#<REMOVE_END>
###############################
echo "Zero field...BSS "
MENU=bss.zerofield;
FILTER=filter_zerofield;
export MENU;
export FILTER;
../run.sh -incmo -utcmo
###############################
echo "Field strength 0.001 ... BSS"
MENU=bss.scf.ccsd.1;
export MENU;
../run.sh -incmo -mw 140 
###############################
echo "Field strength -0.001 ... BSS"
MENU=bss.scf.ccsd.m1;
export MENU;
../run.sh -incmo -mw 140
###############################
rm DFCOEF

