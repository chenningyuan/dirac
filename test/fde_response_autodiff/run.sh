#!/bin/sh

bunzip2 -d *bz2
export MOL=ar.dav2z.mol
export MENU=lda+tf_ll_fdekern.inp
../run.sh -put FRZDNS -put EMBPOT -put EMBGRD
export MENU=lda+tf_ll_nokern.inp
../run.sh -put FRZDNS -put EMBPOT -put EMBGRD

exit 0
