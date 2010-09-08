#!/bin/sh
#
# Set up calculation
#
if [ ! -z "$RUNTYPE" ]; then
   case "${RUNTYPE}" in
      "short")
         MENU=Be.d2h;
         MOL=Be.d2h.mol;
         cp filter_be filter;
         ;;
      "medium")
         MENU=F;
         MOL=F.mol;
         cp filter_f filter;
         ;;
   esac
else
   MENU=Be.d2h;
   MOL=Be.d2h.mol;
   cp filter_be filter;
fi
export MENU;
export MOL;
../run.sh

