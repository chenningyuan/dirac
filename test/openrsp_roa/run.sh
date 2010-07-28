#!/bin/bash

export MOL=h2o2.mol

export MENU=roa_ll.inp
../run.sh -copy DALTON.HES

if [ "$RUNTYPE" = "short" ]; then
   exit 0
fi

export MENU=roa_dc.inp
../run.sh -copy DALTON.HES

exit 0
