#!/bin/sh
#
# Set up calculation
#
if [ ! -z "$RUNTYPE" ]; then
   case "${RUNTYPE}" in
      "short")
         MENU=Be.d2h;
         MOL=Be.d2h.mol;
         cp -f result/filter_short result/filter;
         ;;
      "medium")
         MENU=Be;
         MOL=Be.mol;
         cp -f result/filter_medlong result/filter;
         ;;
      "long")
         MENU=Be;
         MOL=Be.mol;
         cp -f result/filter_medlong result/filter;
         ;;
   esac
else
   MENU=Be.d2h;
   MOL=Be.d2h.mol;
   cp -f result/filter_short result/filter;
fi
export MENU;
export MOL;
../run.sh

