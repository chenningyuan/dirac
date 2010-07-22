#!/bin/sh

export MOL=co.mol

export MENU=response_1z_bz.inp
../run.sh -outcmo -get  PAMXVC

for inp in \
    visual_1z_bz.inp \
    visual_1z_bz_only_para.inp \
    visual_divj.inp \
    visual_divs.inp \
    visual_j.inp \
    visual_jpara.inp \
    visual_rotj.inp \
    visual_rotjpara.inp \
    visual_rots.inp \
    visual_s.inp \
   ;do
    export MENU=$inp
    ../run.sh -incmo -copy PAMXVC
done

# clean up
rm -f DFCOEF* PAMXVC* *tgz

exit 0
