#!/bin/bash

export MOL=ch4.mol

export MENU=pv_beta_nr.inp
../run.sh -copy DALTON.HES

export MENU=pv_gamma_nr.inp
../run.sh -copy DALTON.HES

export MENU=pv_gamma_dc.inp
../run.sh -copy DALTON.HES

# remove restart files

rm -f *tgz

exit 0
