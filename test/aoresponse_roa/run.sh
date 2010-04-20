#!/bin/bash

export MOL=h2o2.mol

export MENU=roa_ll.inp
../run.sh -copy DALTON.HES

#export MENU=roa_ll_restart1.inp
#cp restart1 aoresponse.restart
#../run.sh -copy DALTON.HES -copy aoresponse.restart
#rm -f aoresponse.restart
#
#export MENU=roa_ll_restart2.inp
#cp restart2 aoresponse.restart
#../run.sh -copy DALTON.HES -copy aoresponse.restart
#rm -f aoresponse.restart

if [ "$RUNTYPE" = "medium" ]; then
   exit 0
fi

export MENU=roa_dc.inp
../run.sh -copy DALTON.HES

exit 0
