#!/bin/bash

if [ "$RUNTYPE" = "regular" ]; then
   export MOL=F.mol
   export MENU=f
   cp result/filter_f result/filter
else
   export MOL=Be.d2h.mol
   export MENU=be.d2h
   cp result/filter_be result/filter
fi

../SCRIPTS/run.sh

exit 0
