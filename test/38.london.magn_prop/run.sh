#!/bin/sh
#
# Set up calculation
#
MOL=LiH.lsym.dir;
export MOL;
###############################
echo "levy-leblond NMR shieldings with LAO"
MENU=LiH.levyle.shield_lao_go=cm;
export MENU;
../run.sh
###############################
echo "bss NMR shieldings with LAO"
MENU=LiH.bss.shield_lao_go=cm;
export MENU;
../run.sh
###############################
echo "dc_rkb NMR shieldings with LAO"
MENU=LiH.dc_rkb.shield_lao_go=cm;
export MENU;
../run.sh
