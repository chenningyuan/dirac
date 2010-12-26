#!/bin/sh
#
# Set up calculation
#
if [ ! -z "$RUNTYPE" ]; then
   case "${RUNTYPE}" in
      "short"|"regular")
         MOL=Ne.C2.mol
         ;;
      *)
         MOL=CO.aug-cc-pVDZ.mol
         exit
         ;;
   esac
else
   MOL=CO.aug-cc-pVDZ.mol
fi
MENU=gaunt;
export MENU;
export MOL;
../SCRIPTS/run.sh
