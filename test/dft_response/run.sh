#!/bin/sh

export MOL=hf
export MENU=b3lyp.inp
../run.sh
export MENU=blyp.inp
../run.sh
export MENU=lda.inp
../run.sh

export MOL=he
export MENU=lda_default_dc.inp
../run.sh
export MENU=lda_default_ll.inp
../run.sh
export MENU=lda_default_nr.inp
../run.sh
export MENU=lda_default_x2c.inp
../run.sh
export MENU=lda_default_x2c4.inp
../run.sh

export MOL=he
export MENU=blyp_sdft.inp
../run.sh
export MENU=lda_sdft.inp
../run.sh
export MENU=lda_sdft_collinear.inp
../run.sh
export MENU=lda_nosdft.inp
../run.sh

if [ "$RUNTYPE" = "short" ]; then
   exit
fi

export MOL=hcl
export MENU=blyp_qr_sdft.inp
../run.sh
export MENU=blyp_qr_nosdft.inp
../run.sh
export MENU=lda_qr_sdft.inp
../run.sh

exit 0
