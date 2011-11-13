#!/bin/sh
#
# Set up calculation
#
if [ ! -z "$RUNTYPE" ]; then
   case "${RUNTYPE}" in
      "short")
         MOL=N2.ccpVDZ.mol
         ;;
      "medium")
         MOL=N2.ccpVDZ.mol
         ;;
      "long"|"small_benchmark")
         MOL=N2.ccpVQZ.mol
         ;;
      "benchmark"|"large_benchmark")
         MOL=N2.ccpV5Z.mol
         ;;
   esac
else
   MOL=N2.ccpVTZ.mol
fi
MENU=cc;
export MENU;
export MOL;
../run.sh -mw 130
