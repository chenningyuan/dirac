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
../run.sh -outcmo -get  PAMXVC
export MENU=visual_1z_bz.inp
../run.sh -incmo  -copy PAMXVC
export MENU=visual_1z_bz_only_para.inp
../run.sh -incmo  -copy PAMXVC

export MENU=response_1z_2z.inp
../run.sh -outcmo -get  PAMXVC
export MENU=visual_1z_2z.inp
../run.sh -incmo  -copy PAMXVC

# -----------------------
# visualization with LAOs
# -----------------------

export MOL=co_small.mol
export MENU=shielding_london.inp
#../run.sh -outcmo -get PAMXVC -get AOPROPER
export MENU=visual_shielding_zz_london_natural.inp
#../run.sh -incmo  -copy PAMXVC -copy AOPROPER -get PLOT_3D
#mv PLOT_3D ${MENU}_co.3d
export MENU=visual_shielding_zz_london_natural_onlydirect.inp
#../run.sh -incmo  -copy PAMXVC -copy AOPROPER -get PLOT_3D
#mv PLOT_3D ${MENU}_co.3d
export MENU=visual_shielding_zz_london_natural_onlyresponse.inp
#../run.sh -incmo  -copy PAMXVC -copy AOPROPER -get PLOT_3D
#mv PLOT_3D ${MENU}_co.3d
export MENU=visual_shielding_zz_london_natural_onlyreorth.inp
#../run.sh -incmo  -copy PAMXVC -copy AOPROPER -get PLOT_3D
#mv PLOT_3D ${MENU}_co.3d

# test symmetry
# -------------
for molecule in ne*.mol; do
    export MOL=${molecule%\.mol}

    export MENU=shielding_london.inp
#    ../run.sh -outcmo -get PAMXVC -get AOPROPER

    export MENU=visual_shielding_zz_london_natural.inp
#    ../run.sh -incmo  -copy PAMXVC -copy AOPROPER -get PLOT_3D
#    mv PLOT_3D ${molecule%\.mol}_${MENU}.3d
done


rm -f PAMXVC AOPROPER DFCOEF

export MOL=co_small.mol
export MENU=shielding_london_symcon.inp
#../run.sh -outcmo -get PAMXVC -get AOPROPER
export MENU=visual_shielding_zz_london_symmetric.inp
#../run.sh -incmo  -copy PAMXVC -copy AOPROPER -get PLOT_3D
#mv PLOT_3D ${MENU}_co.3d

export MENU=visual_shielding_zz_london_none.inp
#../run.sh -incmo  -copy PAMXVC -copy AOPROPER -get PLOT_3D
#mv PLOT_3D ${MENU}_co.3d

# test symmetry
# -------------
for molecule in ne*.mol; do
    export MOL=${molecule%\.mol}

    export MENU=shielding_london_symcon.inp
#    ../run.sh -outcmo -get PAMXVC -get AOPROPER

    export MENU=visual_shielding_zz_london_symmetric.inp
#    ../run.sh -incmo  -copy PAMXVC -copy AOPROPER -get PLOT_3D
#    mv PLOT_3D ${molecule%\.mol}_${MENU}.3d

    export MENU=visual_shielding_zz_london_none.inp
#    ../run.sh -incmo  -copy PAMXVC -copy AOPROPER -get PLOT_3D
#    mv PLOT_3D ${molecule%\.mol}_${MENU}.3d
done



### TEST

export MOL=co_small.mol
export MENU=shielding_london.inp
 ../run.sh -outcmo -get PAMXVC -get AOPROPER
export MENU=visual_shielding_zz_london_natural.inp
# ../run.sh -incmo  -copy PAMXVC -copy AOPROPER
export MENU=visual_shielding_zz_london_natural_onlydirect.inp
 ../run.sh -incmo  -copy PAMXVC -copy AOPROPER
export MENU=visual_shielding_zz_london_natural_onlyresponse.inp
# ../run.sh -incmo  -copy PAMXVC -copy AOPROPER


export MOL=co_small_sym.mol
export MENU=shielding_london.inp
 ../run.sh -outcmo -get PAMXVC -get AOPROPER
export MENU=visual_shielding_zz_london_natural.inp
# ../run.sh -incmo  -copy PAMXVC -copy AOPROPER
export MENU=visual_shielding_zz_london_natural_onlydirect.inp
 ../run.sh -incmo  -copy PAMXVC -copy AOPROPER
export MENU=visual_shielding_zz_london_natural_onlyresponse.inp
# ../run.sh -incmo  -copy PAMXVC -copy AOPROPER



export MOL=co_small.mol
export MENU=shielding_london_symcon.inp
# ../run.sh -outcmo -get PAMXVC -get AOPROPER
export MENU=visual_shielding_zz_london_symmetric.inp
# ../run.sh -incmo  -copy PAMXVC -copy AOPROPER
export MENU=visual_shielding_zz_london_symmetric_onlydirect.inp
# ../run.sh -incmo  -copy PAMXVC -copy AOPROPER
export MENU=visual_shielding_zz_london_symmetric_onlyresponse.inp
# ../run.sh -incmo  -copy PAMXVC -copy AOPROPER


export MOL=co_small_sym.mol
export MENU=shielding_london_symcon.inp
# ../run.sh -outcmo -get PAMXVC -get AOPROPER
export MENU=visual_shielding_zz_london_symmetric.inp
# ../run.sh -incmo  -copy PAMXVC -copy AOPROPER
export MENU=visual_shielding_zz_london_symmetric_onlydirect.inp
# ../run.sh -incmo  -copy PAMXVC -copy AOPROPER
export MENU=visual_shielding_zz_london_symmetric_onlyresponse.inp
# ../run.sh -incmo  -copy PAMXVC -copy AOPROPER

export MOL=ne_c2v.mol
export MENU=shielding_london.inp
# ../run.sh -outcmo -get PAMXVC -get AOPROPER
export MENU=visual_shielding_zz_london_natural_onlydirect.inp
# ../run.sh -incmo  -copy PAMXVC -copy AOPROPER
export MENU=visual_shielding_zz_london_natural_onlyresponse.inp
# ../run.sh -incmo  -copy PAMXVC -copy AOPROPER


#export MOL=ne_c1.mol
#export MENU=shielding_london_symcon.inp
# ../run.sh -outcmo -get PAMXVC -get AOPROPER
#export MENU=visual_shielding_zz_london_symmetric_onlydirect.inp
# ../run.sh -incmo  -copy PAMXVC -copy AOPROPER -get PLOT_3D
# mv PLOT_3D ${molecule%\.mol}_${MENU}.3d






exit 0
