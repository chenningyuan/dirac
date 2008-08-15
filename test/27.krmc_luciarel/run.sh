#!/bin/sh
#
# Set up calculation
#
if [ ! -z "$RUNTYPE" ]; then
   case "${RUNTYPE}" in
      "short")
         MENU=Be.d2h;
         MOL=Be.d2h;
         ;;
      "medium")
         MENU=Be;
         MOL=Be;
         ;;
      "long")
         MENU=Be;
         MOL=Be;
         ;;
   esac
else
   MENU=Be.d2h;
   MOL=Be.d2h;
fi
export MENU;
export MOL;
../run.sh
