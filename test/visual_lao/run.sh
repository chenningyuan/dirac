#!/bin/sh


# test symmetry
# -------------

for molecule in co2_*.mol; do
    export MOL=${molecule%\.mol}

    export MENU=response_london_bz.inp
    ../run.sh -outcmo -get PAMXVC -get AOPROPER -get property_info
    export MENU=vis_shield_natural_zz.inp
    ../run.sh -incmo  -copy PAMXVC -copy AOPROPER -copy property_info
    rm -f PAMXVC AOPROPER property_info

    export MENU=response_london_symcon_bz.inp
    ../run.sh -outcmo -get PAMXVC -get AOPROPER -get property_info
    export MENU=vis_shield_symmetric_zz.inp
    ../run.sh -incmo  -copy PAMXVC -copy AOPROPER -copy property_info
    export MENU=vis_shield_none_zz.inp
    ../run.sh -incmo  -copy PAMXVC -copy AOPROPER -copy property_info

done


exit 0
