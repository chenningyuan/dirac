#!/bin/sh
#
# Set up calculation
#
if [ ! -z "$RUNTYPE" ]; then
   case "${RUNTYPE}" in
      "short")
         MOL=N2.ccpVDZ
         ;;
      "medium")
         MOL=N2.ccpVDZ
         ;;
      "long"|"small_benchmark")
         MOL=N2.ccpVQZ
         ;;
      "benchmark"|"large_benchmark")
         MOL=N2.ccpV5Z
         ;;
   esac
else
   MOL=N2.ccpVTZ
fi
MENU=cc;
export MENU;
export MOL;
../run.sh
