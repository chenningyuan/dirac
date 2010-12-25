#!/bin/bash

export MENU=acetone
export MOL=acetone
../SCRIPTS/run.sh

export MENU=HeH+
export MOL=HeH+
../SCRIPTS/run.sh

if [ "$RUNTYPE" = "medium" ]; then
   export MENU=HCCF.x2c_zdipf.scf_mp2_geopt.inp
   export MOL=HCCF.sto-2g.mol
   ../SCRIPTS/run.sh

   export MENU=1.3.5_trisilacyclohexane
   export MOL=1.3.5_trisilacyclohexane.mol
   ../SCRIPTS/run.sh
fi    

exit 0
