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
../run.sh -utcmo
###############################
echo "Field strength 0.001 ... DC"
MENU=scf.ccsd.1;
export MENU;
../run.sh -incmo
###############################
echo "Field strength -0.001 ... DC"
MENU=scf.ccsd.m1;
export MENU;
../run.sh -incmo
###############################
rm DFCOEF
###############################
###############################
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
../run.sh -incmo
###############################
echo "Field strength -0.001 ... BSS"
MENU=bss.scf.ccsd.m1;
export MENU;
../run.sh -incmo
###############################
rm DFCOEF

