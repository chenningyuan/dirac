#!/bin/bash

for mol in ar.av2z.mol ar.av3z.mol ar.av4z.mol ar.av5z.mol ar.v2z.mol ar.v3z.mol ar.v4z.mol ar.v5z.mol ; do
   for inp in lda_nonrel_fdekern.inp lda_nonrel_nofde.inp lda_nonrel_nokern.inp ; do

      ../../../pam -mpi 1 -put EMBPOT -put FRZDN1 -put EMBGRD -clean $mol $inp 
   done
done
