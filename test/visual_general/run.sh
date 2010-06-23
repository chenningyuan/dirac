#!/bin/sh

#for molecule in ne*mol; do
for molecule in ne_auto.mol; do
    export MOL=${molecule%\.mol}

    export MENU=response_ez_ez.inp
#    ../run.sh -outcmo -get  PAMXVC
    export MENU=visual_zo.inp
#    ../run.sh -incmo  -copy PAMXVC -get PLOT_3D
#    mv PLOT_3D ${molecule%\.mol}_zo.3d
    export MENU=visual_ez_ez.inp
#    ../run.sh -incmo  -copy PAMXVC

    export MENU=response_bz_bz.inp
#    ../run.sh -outcmo -get  PAMXVC
    export MENU=visual_bz_bz.inp
#    ../run.sh -incmo  -copy PAMXVC
done

export MOL=co.mol

export MENU=response_1z_bz.inp
#../run.sh -outcmo -get  PAMXVC
export MENU=visual_1z_bz.inp
#../run.sh -incmo  -copy PAMXVC
export MENU=visual_1z_bz_only_para.inp
#../run.sh -incmo  -copy PAMXVC

export MENU=response_1z_2z.inp
#../run.sh -outcmo -get  PAMXVC
export MENU=visual_1z_2z.inp
#../run.sh -incmo  -copy PAMXVC


export MENU=shielding_london.inp
../run.sh -outcmo -get PAMXVC -get AOPROPER
export MENU=visual_shielding_xx_london.inp
../run.sh -incmo  -copy PAMXVC -copy AOPROPER -get PLOT_3D
mv PLOT_3D ${MENU}_co.3d
export MENU=visual_shielding_yy_london.inp
#../run.sh -incmo  -copy PAMXVC -copy AOPROPER -get PLOT_3D
#mv PLOT_3D ${MENU}_co.3d
export MENU=visual_shielding_zz_london.inp
#../run.sh -incmo  -copy PAMXVC -copy AOPROPER -get PLOT_3D
#mv PLOT_3D ${MENU}_co.3d
export MENU=visual_current_bx_london.inp
#../run.sh -incmo  -copy PAMXVC -copy AOPROPER -get PLOT_2D
#mv PLOT_2D ${MENU}_co.2d

export MENU=shielding_london_symcon.inp
#../run.sh -outcmo -get PAMXVC -get AOPROPER
export MENU=visual_shielding_xx_london_symcon.inp
#../run.sh -incmo  -copy PAMXVC -copy AOPROPER -get PLOT_3D
#mv PLOT_3D ${MENU}_co.3d
export MENU=visual_shielding_yy_london.inp
#../run.sh -incmo  -copy PAMXVC -copy AOPROPER -get PLOT_3D
#mv PLOT_3D ${MENU}_co.3d
export MENU=visual_shielding_zz_london.inp
#../run.sh -incmo  -copy PAMXVC -copy AOPROPER -get PLOT_3D
#mv PLOT_3D ${MENU}_co.3d
export MENU=visual_current_bx_london.inp
#../run.sh -incmo  -copy PAMXVC -copy AOPROPER -get PLOT_2D
#mv PLOT_2D ${MENU}_co.2d

export MENU=visual_shielding_xx_london_none.inp
#../run.sh -incmo  -copy PAMXVC -copy AOPROPER -get PLOT_3D
#mv PLOT_3D ${MENU}_co.3d

exit 0
