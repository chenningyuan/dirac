#!/bin/sh
#
# Set up calculation
#
MOL=LiH.lsym.dir;
FILTER=filter;
export FILTER;
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
###############################
echo "levy-leblond magnetic susceptibilities with LAO"
FILTER=filter_susc;
export FILTER;
MENU=LiH.levyle.suscp_lao_go=cm;
export MENU;
../run.sh
###############################

