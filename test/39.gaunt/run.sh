#!/bin/sh
#
# Set up calculation
#
if [ ! -z "$RUNTYPE" ]; then
   case "${RUNTYPE}" in
      "short"|"medium")
         MOL=Ne.C2
         ;;
      *)
         MOL=CO.aug-cc-pVTZ
         exit
         ;;
   esac
else
   MOL=CO.aug-cc-pVTZ
fi
MENU=gaunt;
export MENU;
export MOL;
../run.sh
