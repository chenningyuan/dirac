#!/bin/sh
#
# Set up calculation
#
if [ ! -z "$RUNTYPE" ]; then
   case "${RUNTYPE}" in
      "short")
         MENU=Be.gascip.bss+mfsso.krmc;
         MOL=Be.gascip.mol;
         ;;
      "medium")
         MENU=Be;
         MOL=Be.mol;
         ;;
      "long")
         MENU=Be.gascip;
         MOL=Be.gascip.mol;
         ;;
   esac
else
   MENU=Be.gascip.bss+mfsso.krmc;
   MOL=Be.gascip.mol;
fi
export MENU;
export MOL;
../run.sh

