#!/bin/sh

export MOL=co.mol

export MENU=response_1z_bz.inp
../run.sh -outcmo -get  PAMXVC

for inp in visual*inp; do
    export MENU=$inp
    ../run.sh -incmo -copy PAMXVC
done

# clean up
rm -f DFCOEF* PAMXVC* *tgz

exit 0
