#!/bin/sh

export MOL=hf

export MENU=b3lyp.inp
../SCRIPTS/run.sh

export MENU=blyp.inp
../SCRIPTS/run.sh

export MENU=lda.inp
../SCRIPTS/run.sh

export MENU=pbex.inp
../SCRIPTS/run.sh

export MENU=pbec.inp
../SCRIPTS/run.sh

export MENU=pbe0.inp
../SCRIPTS/run.sh

export MENU=pbe.inp
../SCRIPTS/run.sh

export MOL=he

export MENU=blyp_sdft.inp
../SCRIPTS/run.sh

export MENU=lda_sdft.inp
../SCRIPTS/run.sh

exit 0
