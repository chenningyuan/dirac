#!/bin/sh

export MOL=fh

export MENU=lda.inp
../run.sh -dftmem 100

export MENU=blyp.inp
../run.sh -dftmem 100

exit 0
