#!/bin/sh

bunzip2 -d *bz2
export MOL=water.sto6g.mol
export MENU=fde_ll_lda.inp
../run.sh -put h2o_in_2h2o.frzdens -put h2o_in_2h2o.vemb
export MENU=nofde_ll_lda.inp
../run.sh

export MENU=fde_ll_bp86.inp
../run.sh -put h2o_in_2h2o.frzdens -put h2o_in_2h2o.vemb
export MENU=nofde_ll_bp86.inp
../run.sh

exit 0
