#!/bin/bash

export MOL=h2o2.mol

export MENU=nr.inp
../SCRIPTS/run.sh

export MOL=hf.mol

# test nonrel transformation

export MENU=nr.inp
../SCRIPTS/run.sh
export MENU=nr_shifted.inp
../SCRIPTS/run.sh

# test nonrel limit

export MENU=dc_skipep.inp
../SCRIPTS/run.sh
export MENU=dc_skipep_cval=800.inp
../SCRIPTS/run.sh

# test rel transformation

export MENU=dc.inp
../SCRIPTS/run.sh
export MENU=dc_shifted.inp
../SCRIPTS/run.sh

# remove restart files

rm -f *tgz

exit 0
