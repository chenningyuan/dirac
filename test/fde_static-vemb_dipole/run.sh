#!/bin/sh

bzip2 -d EMBPOT.k.nr.acc4.lda+tfk.fdem_h2o+2h2o_tz2p_nosym.bz2

ln -sf EMBPOT.k.nr.acc4.lda+tfk.fdem_h2o+2h2o_tz2p_nosym vemb.grid

for mol in water.sto6g water.sto6g.c2v; do
   export MOL=$mol.mol
   for inp in fde_dc_lda fde_ll_lda fde_twocomp_lda fde_dc_cc; do
      export MENU=$inp.inp
      ../run.sh -mw 80 -nw 80 -put vemb.grid 
   done
done

bzip2 -9 EMBPOT.k.nr.acc4.lda+tfk.fdem_h2o+2h2o_tz2p_nosym

exit 0

