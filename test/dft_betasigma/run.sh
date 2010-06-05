#!/bin/sh

export MOL=he
export MENU=lda_sigma.inp
../run.sh
export MENU=lda_betasigma.inp
../run.sh
export MENU=blyp_sigma.inp
../run.sh
export MENU=blyp_betasigma.inp
../run.sh

exit 0
