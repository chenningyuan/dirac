#!/bin/sh

for molecule in ne*mol; do
    export MOL=${molecule%\.mol}

    export MENU=response_ez_ez.inp
    ../run.sh -outcmo -get  PAMXVC
    export MENU=visual_zo.inp
    ../run.sh -incmo  -copy PAMXVC -get PLOT_3D
    mv PLOT_3D ${molecule%\.mol}_zo.3d
    export MENU=visual_ez_ez.inp
    ../run.sh -incmo  -copy PAMXVC

    export MENU=response_bz_bz.inp
    ../run.sh -outcmo -get  PAMXVC
    export MENU=visual_bz_bz.inp
    ../run.sh -incmo  -copy PAMXVC
done

export MOL=co.mol

export MENU=response_1z_bz.inp
../run.sh -outcmo -get  PAMXVC
export MENU=visual_1z_bz.inp
../run.sh -incmo  -copy PAMXVC
export MENU=visual_1z_bz_only_para.inp
../run.sh -incmo  -copy PAMXVC

export MENU=response_1z_2z.inp
../run.sh -outcmo -get  PAMXVC
export MENU=visual_1z_2z.inp
../run.sh -incmo  -copy PAMXVC

# clean up
rm -f DFCOEF* PAMXVC*

exit 0
