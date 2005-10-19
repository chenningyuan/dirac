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
         MOL=CO.aug-cc-pVDZ
         exit
         ;;
   esac
else
   MOL=CO.aug-cc-pVDZ
fi
MENU=gaunt;
export MENU;
export MOL;
../run.sh
