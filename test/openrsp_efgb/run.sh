#!/bin/bash

export MOL=h2o2.mol

export MENU=nr_hf.inp
../run.sh

export MOL=fh.mol

# test nonrel transformation

export MENU=nr_hf.inp
../run.sh
export MENU=nr_hf_shifted.inp
../run.sh

# test nonrel limit

export MENU=dc_hf_skipep.inp
../run.sh
export MENU=dc_hf_skipep_cval=3000.inp
../run.sh

# test rel transformation

export MENU=dc_hf.inp
../run.sh
export MENU=dc_hf_shifted.inp
../run.sh

# remove restart files

rm -f *tgz

exit 0
