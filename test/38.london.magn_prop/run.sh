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


MOL=LiH.lsym.dir.mol;
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
###############################
if [ "$RUNTYPE" = "short" ]; then
   echo
   echo "Completed short test"
   echo "<-----"
   exit
fi
###############################
echo "bss+mfsso NMR shieldings with LAO"
MENU=LiH.bss+mfsso.shield_lao_go=cm;
export MENU;
../run.sh
###############################
echo "bss_sfB NMR shieldings with LAO"
MENU=LiH.bss_sfB.shield_lao_go=cm;
export MENU;
../run.sh
###############################
echo "dc_rkb NMR shieldings with LAO"
MENU=LiH.dc_rkb.shield_lao_go=cm;
export MENU;
../run.sh
###############################
echo "dc_sf NMR shieldings with LAO"
MENU=LiH.dc_sf.shield_lao_go=cm;
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

