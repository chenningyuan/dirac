#!/bin/sh
#
# Set up calculation
#
if [ ! -z "$RUNTYPE" ]; then
   case "${RUNTYPE}" in
      "short")
         MOL=H2O
         ;;
      "medium")
         MOL=CH4
         ;;
      "long"|"small_benchmark")
         MOL=SF6
         ;;
      "benchmark"|"large_benchmark")
         MOL=SF6
         ;;
   esac
else
   MOL=H2O
fi
MENU=b3lyp;
export MENU;
export MOL;
../run.sh
