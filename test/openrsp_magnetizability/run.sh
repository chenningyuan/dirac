#!/bin/bash

export MOL=fh.mol

# test nonrel transformation

export MENU=nr_hf.inp
../SCRIPTS/run.sh
export MENU=nr_hf_shifted.inp
../SCRIPTS/run.sh

# test nonrel limit

export MENU=dc_hf_skipep.inp
../SCRIPTS/run.sh
export MENU=dc_hf_skipep_cval=900.inp
../SCRIPTS/run.sh

# test rel transformation

export MENU=dc_hf_cval=30.inp
../SCRIPTS/run.sh
export MENU=dc_hf_shifted_cval=30.inp
../SCRIPTS/run.sh

# remove restart files

rm -f *tgz

exit 0
