#!/bin/sh

export MOL=hehe
export MENU=lda_levyle_nosdft_london_shield.inp
../run.sh
export MENU=blyp_levyle_nosdft_london_shield.inp
../run.sh
exit 0
export MENU=b3lyp_levyle_nosdft_london_shield.inp
../run.sh


export MOL=hehe_sym
export MENU=lda_levyle_nosdft_london_shield.inp
../run.sh
export MENU=blyp_levyle_nosdft_london_shield.inp
../run.sh
export MENU=b3lyp_levyle_nosdft_london_shield.inp
../run.sh

exit 0
