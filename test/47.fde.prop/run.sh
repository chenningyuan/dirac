#!/bin/sh

bunzip2 -d *bz2

export MOL=water.sto6g.mol

export MENU=fde_ll_lda.inp
../run.sh -put vemb.grid
export MENU=fde_dc_lda.inp
../run.sh -put vemb.grid
export MENU=fde_twocomp_lda.inp
../run.sh -put vemb.grid
export MENU=fde_ll_cc.inp
../run.sh -put vemb.grid
export MENU=fde_dc_cc.inp
../run.sh -put vemb.grid

if [ ! -z "$RUNTYPE" ]; then
   case "${RUNTYPE}" in
      "medium"|"long")
      export MENU=fde_twocomp_cc.inp
       ../run.sh -mw 80 -nw 80 -put vemb.grid
         ;;
   esac
fi



export MOL=water.sto6g.c2v.mol

export MENU=fde_ll_lda.inp
../run.sh -put vemb.grid
export MENU=fde_dc_lda.inp
../run.sh -put vemb.grid
export MENU=fde_twocomp_lda.inp
../run.sh -put vemb.grid
export MENU=fde_ll_cc.inp
../run.sh -put vemb.grid
export MENU=fde_dc_cc.inp
../run.sh -put vemb.grid


if [ ! -z "$RUNTYPE" ]; then
   case "${RUNTYPE}" in
      "medium")
      export MENU=fde_twocomp_cc.inp
       ../run.sh -mw 80 -nw 80 -put vemb.grid
         ;;
   esac
fi

#export MOL=water.sto6g.mol
#
#export MENU=nofde_ll_lda.inp
#../run.sh 
#export MENU=nofde_dc_lda.inp
#../run.sh 
#export MENU=nofde_twocomp_lda.inp
#../run.sh 
#
#
#export MOL=water.sto6g.c2v.mol
#
#export MENU=nofde_ll_lda.inp
#../run.sh 
#export MENU=nofde_dc_lda.inp
#../run.sh 
#export MENU=nofde_twocomp_lda.inp
#../run.sh 


exit 0

