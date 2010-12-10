#!/bin/sh

export MOL=he
export MENU=lda_sigma.inp
../SCRIPTS/run.sh
export MENU=lda_betasigma.inp
../SCRIPTS/run.sh
export MENU=blyp_sigma.inp
../SCRIPTS/run.sh
export MENU=blyp_betasigma.inp
../SCRIPTS/run.sh

exit 0
