#!/bin/sh

export MENU=lda_levyle_nosdft_london_shield.inp
export MOL=hehe
../run.sh
export MOL=hehe_sym
../run.sh

exit 0
