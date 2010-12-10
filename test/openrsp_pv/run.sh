#!/bin/bash

export MOL=ch4.mol

export MENU=pv_beta_nr.inp
../SCRIPTS/run.sh -copy DALTON.HES

export MENU=pv_gamma_nr.inp
../SCRIPTS/run.sh -copy DALTON.HES

export MENU=pv_gamma_dc.inp
../SCRIPTS/run.sh -copy DALTON.HES

# remove restart files

rm -f *tgz

exit 0
