#!/bin/sh
#
# Set up calculation
#
if [ ! -z "$RUNTYPE" ]; then
   case "${RUNTYPE}" in
      "short")
         export MENU=HeH+;
         export MOL=HeH+;
         ../SCRIPTS/run.sh
         ;;
      "medium")
         export MENU=acetone;
         export MOL=acetone;
         ../SCRIPTS/run.sh
         export MENU=HCCF.x2c_zdipf.scf_mp2_geopt.inp;
         export MOL=HCCF.sto-2g.mol;
          ../SCRIPTS/run.sh
         ;;
      "long")
         export MENU=1.3.5_trisilacyclohexane;
         export MOL=1.3.5_trisilacyclohexane.mol;
         ../SCRIPTS/run.sh
         ;;
   esac
else
# if no runtype test, use acetone as the default
   export MENU=acetone;
   export MOL=acetone;
   ../SCRIPTS/run.sh
fi
