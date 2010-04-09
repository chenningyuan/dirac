#!/bin/sh

export MOL=hf
export MENU=blyp.inp
../run.sh -dftmem 10
export MENU=lda.inp
../run.sh -dftmem 10

exit 0
