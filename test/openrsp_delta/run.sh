#!/bin/sh

export MOL=fh

export MENU=lda.inp
../SCRIPTS/run.sh

export MENU=blyp.inp
../SCRIPTS/run.sh

exit 0
