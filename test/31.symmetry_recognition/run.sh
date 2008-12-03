#!/bin/sh

export MENU=lda

export MOL=H2O
../run.sh

export MOL=CH4
../run.sh

if [ "$RUNTYPE" = "short" ]; then
   echo
   echo "Completed short test"
   echo "<-----"
   exit
fi

export MOL=SF6
../run.sh

exit 0
