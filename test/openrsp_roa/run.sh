#!/bin/bash

export MOL=h2o2.mol

export MENU=roa_ll.inp
../SCRIPTS/run.sh --copy=DALTON.HES

if [ "$RUNTYPE" = "short" ]; then
   exit 0
fi

export MENU=roa_dc.inp
../SCRIPTS/run.sh --copy=DALTON.HES

exit 0
