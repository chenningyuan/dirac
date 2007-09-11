#!/bin/sh
#
# Set up calculation
#
if [ ! -z "$RUNTYPE" ]; then
   case "${RUNTYPE}" in
      "short")
         MENU=ecp_d2h
         ;;
      "medium")
         MENU=ecp_d2h
         ;;
      "long")
         MENU=ecp_d2h
         ;;
   esac
else
   MENU=ecp_d2h;
fi
export MENU;
MOL=f_arpp_d2h; export MOL;
../run.sh
MOL=f_sopp_d2h; export MOL;
../run.sh


