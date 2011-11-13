#!/bin/sh
#
# Set up calculation
#
MOL=HF.p-ccpvdz;
export MOL;
###############################
echo "Zero field..."
MENU=zerofield;
FILTER=filter_zerofield;
export MENU;
export FILTER;
../run.sh -utcmo
###############################
echo "Field strength 0.001"
MENU=scf.ccsd.1;
export MENU;
../run.sh -incmo
###############################
echo "Field strength -0.001"
MENU=scf.ccsd.m1;
export MENU;
../run.sh -incmo
###############################
rm DFCOEF
