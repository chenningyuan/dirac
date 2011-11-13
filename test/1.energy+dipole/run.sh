#!/bin/sh
#
# Set up calculation
#
if [ ! -z "$RUNTYPE" ]; then
   case "${RUNTYPE}" in
      "short")
         MENU=ccsd.small
         ;;
      "medium")
         MENU=ccsd
         ;;
      "long")
         MENU=ccsd
         ;;
   esac
else
   MENU=ccsd;
fi
MOL=H2O;
export MENU;
export MOL;
../run.sh
