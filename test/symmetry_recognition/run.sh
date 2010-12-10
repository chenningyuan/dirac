#!/bin/sh

export MENU=lda

export MOL=H2O
../SCRIPTS/run.sh

export MOL=CH4
../SCRIPTS/run.sh

if [ "$RUNTYPE" = "short" ]; then
   echo
   echo "Completed short test"
   echo "<-----"
   exit
fi

export MOL=SF6
../SCRIPTS/run.sh

exit 0
