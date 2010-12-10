#!/bin/sh

bzip2 -d ../fde_static-vemb_dipole/EMBPOT.k.nr.acc4.lda+tfk.fdem_h2o+2h2o_tz2p_nosym.bz2

ln -sf ../fde_static-vemb_dipole/EMBPOT.k.nr.acc4.lda+tfk.fdem_h2o+2h2o_tz2p_nosym vemb.grid
ln -sf ../fde_static-vemb_dipole/EMBPOT.k.nr.acc4.lda+tfk.fdem_h2o+2h2o_tz2p_nosym GRIDOUT 

for mol in water.sto6g water.sto6g.c2v; do
   export MOL=$mol.mol
   for inp in fde_dc_hf_exportdens fde_dc_lda_exportdens fde_dc_mp2_exportdens; do
      export MENU=$inp.inp
      ../SCRIPTS/run.sh -mw 80 -nw 80 -put vemb.grid -put GRIDOUT 
   done
done

bzip2 -9 ../fde_static-vemb_dipole/EMBPOT.k.nr.acc4.lda+tfk.fdem_h2o+2h2o_tz2p_nosym

exit 0
