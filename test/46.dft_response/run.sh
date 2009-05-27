#!/bin/sh

export MOL=hf
export MENU=b3lyp.inp
../run.sh
export MENU=blyp.inp
../run.sh
export MENU=lda.inp
../run.sh

export MOL=he
export MENU=blyp_sdft.inp
../run.sh
export MENU=lda_sdft.inp
../run.sh

export MOL=hcl
export MENU=blyp_qr_sdft.inp
#../run.sh
export MENU=lda_qr_sdft.inp
#../run.sh

exit 0
