#!/bin/bash

export MOL=mol.mol

export MENU=lda.inp
../SCRIPTS/run.sh
export MENU=blyp.inp
../SCRIPTS/run.sh
export MENU=blyp_shifted.inp
../SCRIPTS/run.sh

rm -f *tgz

exit 0
