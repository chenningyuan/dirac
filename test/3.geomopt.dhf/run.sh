#!/bin/sh
#
# Set up calculation
#
if [ ! -z "$RUNTYPE" ]; then
   case "${RUNTYPE}" in
      "short")
         MENU=HeH+;
         MOL=HeH+;
         ;;
      "medium")
         MENU=acetone;
         MOL=acetone;
         ;;
      "long")
         MENU=1.3.5_trisilacyclohexane;
         MOL=1.3.5_trisilacyclohexane;
         ;;
   esac
else
   MENU=acetone;
   MOL=acetone;
fi
export MENU;
export MOL;
../run.sh
