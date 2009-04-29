#!/bin/sh

export MOL=hf

export MENU=b3lyp.inp
../run.sh

export MENU=blyp.inp
../run.sh

export MENU=lda.inp
../run.sh

exit 0
