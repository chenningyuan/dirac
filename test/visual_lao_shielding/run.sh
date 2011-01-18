#!/bin/sh


export MOL=h2o


# Levy-Leblond HF
# ---------------
export MENU=ll_symmetric
../SCRIPTS/run.sh --outcmo --get="PAMXVC AOPROPER property_info"
export MENU=vis_sigmazz_ll_symmetric
../SCRIPTS/run.sh --incmo  --copy="PAMXVC AOPROPER property_info"
export MENU=vis_sigmazz_ll_none
../SCRIPTS/run.sh --incmo  --copy="PAMXVC AOPROPER property_info"
rm -f DFCOEF PAMXVC AOPROPER property_info

export MENU=ll_natural
../SCRIPTS/run.sh --outcmo --get="PAMXVC AOPROPER property_info"
export MENU=vis_sigmazz_ll_natural
../SCRIPTS/run.sh --incmo  --copy="PAMXVC AOPROPER property_info"
rm -f DFCOEF PAMXVC AOPROPER property_info


# 4c DHF
# ------
export MENU=dc_symmetric
../SCRIPTS/run.sh --outcmo --get="PAMXVC AOPROPER property_info"
export MENU=vis_sigmazz_dc_symmetric
../SCRIPTS/run.sh --incmo  --copy="PAMXVC AOPROPER property_info"
export MENU=vis_sigmazz_dc_none
../SCRIPTS/run.sh --incmo  --copy="PAMXVC AOPROPER property_info"
rm -f DFCOEF PAMXVC AOPROPER property_info

export MENU=dc_natural
../SCRIPTS/run.sh --outcmo --get="PAMXVC AOPROPER property_info"
export MENU=vis_sigmazz_dc_natural
../SCRIPTS/run.sh --incmo  --copy="PAMXVC AOPROPER property_info"
rm -f DFCOEF PAMXVC AOPROPER property_info


exit 0
