#!/bin/bash

if [ "$RUNTYPE" = "regular" ]; then
   export MOL=H3
   export MENU=cosci.H3
else
   export MOL=Li2
   export MENU=cosci.Li2
fi

../SCRIPTS/run.sh --mw=120

exit 0
