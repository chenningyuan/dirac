#!/bin/sh

export MOL=hf
export MENU=b3lyp.inp
../SCRIPTS/run.sh
export MENU=blyp.inp
../SCRIPTS/run.sh
export MENU=lda.inp
../SCRIPTS/run.sh

export MOL=he
export MENU=lda_default_dc.inp
../SCRIPTS/run.sh
export MENU=lda_default_ll.inp
../SCRIPTS/run.sh
export MENU=lda_default_nr.inp
../SCRIPTS/run.sh
export MENU=lda_default_x2c.inp
../SCRIPTS/run.sh
export MENU=lda_default_x2c4.inp
../SCRIPTS/run.sh

export MOL=he
export MENU=blyp_sdft.inp
../SCRIPTS/run.sh
export MENU=lda_sdft.inp
../SCRIPTS/run.sh
export MENU=lda_sdft_collinear.inp
../SCRIPTS/run.sh
export MENU=lda_nosdft.inp
../SCRIPTS/run.sh

if [ "$RUNTYPE" = "short" ]; then
   exit
fi

export MOL=hcl
export MENU=blyp_qr_sdft.inp
../SCRIPTS/run.sh
export MENU=blyp_qr_nosdft.inp
../SCRIPTS/run.sh
export MENU=lda_qr_sdft.inp
../SCRIPTS/run.sh

exit 0
