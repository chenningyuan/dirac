#!/bin/bash

export MOL=h2o2.mol

export MENU=nr.inp
../run.sh

export MOL=hf.mol

# test nonrel transformation

export MENU=nr.inp
../run.sh
export MENU=nr_shifted.inp
../run.sh

# test nonrel limit

export MENU=dc_skipep.inp
../run.sh
export MENU=dc_skipep_cval=800.inp
../run.sh

# test rel transformation

export MENU=dc.inp
../run.sh
export MENU=dc_shifted.inp
../run.sh

# remove restart files

rm -f *tgz

exit 0
