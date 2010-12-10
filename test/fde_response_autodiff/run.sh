#!/bin/sh

bzip2 -d FRZDNS.k.nr.acc4.lda+tfk.fdem_h2o+2h2o_tz2p_nosym.bz2
bzip2 -d ../fde_static-vemb_dipole/EMBPOT.k.nr.acc4.lda+tfk.fdem_h2o+2h2o_tz2p_nosym.bz2

ln -sf FRZDNS.k.nr.acc4.lda+tfk.fdem_h2o+2h2o_tz2p_nosym FRZDNS
ln -sf ../fde_static-vemb_dipole/EMBPOT.k.nr.acc4.lda+tfk.fdem_h2o+2h2o_tz2p_nosym EMBPOT 

export MENU=lda+tf_ll_5_fdekern.inp

export MOL=h2o.6311++G_2d2p.c2v.mol
../SCRIPTS/run.sh -put FRZDNS -put EMBPOT 
export MOL=h2o.6311++G_2d2p.mol
../SCRIPTS/run.sh -put FRZDNS -put EMBPOT 

bzip2 -9 FRZDNS.k.nr.acc4.lda+tfk.fdem_h2o+2h2o_tz2p_nosym 
bzip2 -9 ../fde_static-vemb_dipole/EMBPOT.k.nr.acc4.lda+tfk.fdem_h2o+2h2o_tz2p_nosym 

exit 0
