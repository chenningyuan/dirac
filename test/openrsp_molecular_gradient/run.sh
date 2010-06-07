#!/bin/bash

export MOL=h2o2.mol
for inp in *inp; do
    export MENU=$inp
    ../run.sh
done

exit 0

export MOL=c3h6o.mol
for inp in hf_dc.inp hf_dc_cval=5.inp hf_ll.inp; do
    export MENU=$inp
    ../run.sh
done

exit 0
