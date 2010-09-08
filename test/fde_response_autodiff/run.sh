#!/bin/sh


export MENU=lda+tf_ll_fdekern.inp

bunzip2 -d *.k.nr.acc10.lda+tf.fdem_ar+ar_qz4p.bz2

ln -sf FRZDNS.k.nr.acc10.lda+tf.fdem_ar+ar_qz4p FRZDNS
ln -sf EMBPOT.k.nr.acc10.lda+tf.fdem_ar+ar_qz4p EMBPOT

export MOL=ar.dav2z.mol

../run.sh -put FRZDNS -put EMBPOT 

bzip2 -9 *.k.nr.acc10.lda+tf.fdem_ar+ar_qz4p

#
# doing water embedded in the potential of two other water
# molecules
# 

bunzip2 -d *.k.nr.acc4.lda+tf.fdem_h2o+2h2o_dzp.bz2

ln -sf FRZDNS.k.nr.acc4.lda+tf.fdem_h2o+2h2o_dzp FRZDNS
ln -sf EMBPOT.k.nr.acc4.lda+tf.fdem_h2o+2h2o_dzp EMBPOT

export MOL=h2o.dav2z.c2v.mol
../run.sh -put FRZDNS -put EMBPOT 
export MOL=h2o.dav2z.mol
../run.sh -put FRZDNS -put EMBPOT 

bzip2 -9 *.k.nr.acc4.lda+tf.fdem_h2o+2h2o_dzp

exit 0
