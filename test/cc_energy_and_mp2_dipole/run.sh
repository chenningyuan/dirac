#!/bin/bash

export MOL=H2O

if [ "$RUNTYPE" = "regular" ]; then
   export MENU=ccsd
else
   export MENU=ccsd.small
fi

../SCRIPTS/run.sh --mw=120 --aw=130

exit 0
