#!/bin/bash

if [ "$RUNTYPE" = "medium" ]; then
   export MOL=F.mol
   export MENU=F
   cp result/filter_f result/filter
else
   export MOL=Be.d2h.mo
   export MENU=Be.d2h
   cp result/filter_be result/filter
fi

../SCRIPTS/run.sh

exit 0
