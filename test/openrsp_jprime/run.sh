#!/bin/bash

export MOL=hf.mol

export MENU=nr.inp
../run.sh
export MENU=sf_skipep.inp
../run.sh
export MENU=dc_skipep.inp
../run.sh
export MENU=dc_skipep_cval=800.inp
../run.sh

rm -f *tgz

exit 0
