#!/bin/sh
#
# Set up calculation
#
if [ ! -z "$RUNTYPE" ]; then
   case "${RUNTYPE}" in
      "short")
         MENU=cosci.Li2;
	   MOL=Li2;
         ;;
      "medium")
         MENU=cosci.H3;
	   MOL=H3;
         ;;
      "long")
         MENU=cosci.H3;
	   MOL=H3;
         ;;
   esac
else
   MENU=cosci.Li2;
   MOL=Li2;
fi
export MENU;
export MOL;
../SCRIPTS/run.sh -mw 120
