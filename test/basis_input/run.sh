#!/bin/sh
#
# Set up calculation
#  
MENU=ch4; export MENU;
MOL=ch4_temp; export MOL;
../SCRIPTS/run.sh
MOL=ch4_explicit; export MOL;
../SCRIPTS/run.sh
MOL=ch4_family; export MOL;
../SCRIPTS/run.sh
MOL=ch4_dualfamily; export MOL;
../SCRIPTS/run.sh
MOL=ch4_basis; export MOL;
#../SCRIPTS/run.sh --basdir=`pwd`
../SCRIPTS/run.sh 
MOL=ch4_bsse; export MOL;
#../SCRIPTS/run.sh --basdir=`pwd`
../SCRIPTS/run.sh 
MOL=ch4_basis_vtz; export MOL;
#../SCRIPTS/run.sh --basdir=`pwd`
../SCRIPTS/run.sh 
MOL=ch4_multibasis_vtz; export MOL;
#../SCRIPTS/run.sh --basdir=`pwd`
../SCRIPTS/run.sh 
MOL=ch4_molfbas; export MOL;
#../SCRIPTS/run.sh --basdir=`pwd`
../SCRIPTS/run.sh 
