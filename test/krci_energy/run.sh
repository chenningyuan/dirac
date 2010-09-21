#!/bin/sh
#
# Set up calculation
#
if [ ! -z "$RUNTYPE" ]; then
   case "${RUNTYPE}" in
      "short")
         MENU=Be.d2h;
         MOL=Be.d2h.mol;
         cp result/filter_be result/filter;
         ;;
      "medium")
         MENU=F;
         MOL=F.mol;
         cp result/filter_f result/filter;
         ;;
      "long")
         MENU=F;
         MOL=F.mol;
         cp result/filter_f result/filter;
         ;;
   esac
else
   MENU=Be.d2h;
   MOL=Be.d2h.mol;
   cp result/filter_be result/filter;
fi
export MENU;
export MOL;
../run.sh

