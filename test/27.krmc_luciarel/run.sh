#!/bin/sh
#
# Set up calculation
#
if [ ! -z "$RUNTYPE" ]; then
   case "${RUNTYPE}" in
      "short")
         MENU=Be.d2h;
         MOL=Be.d2h.mol;
         ;;
      "medium")
         MENU=Be;
         MOL=Be.mol;
         ;;
      "long")
         MENU=Be;
         MOL=Be.mol;
         ;;
   esac
else
   MENU=Be.d2h;
   MOL=Be.d2h.mol;
fi
export MENU;
export MOL;
../run.sh
