#!/bin/bash

export MOL=mol.mol

export MENU=lda.inp
../run.sh
export MENU=blyp.inp
../run.sh
export MENU=blyp_shifted.inp
../run.sh

rm -f *tgz

exit 0
