#!/bin/sh

export MOL=fh

export MENU=lda.inp
../run.sh

export MENU=blyp.inp
../run.sh

exit 0
