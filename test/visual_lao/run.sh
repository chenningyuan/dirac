#!/bin/sh

if [ "$RUNTYPE" = "short" ]; then
   for molecule in co2_d2h.mol; do
       export MOL=${molecule%\.mol}
   
       export MENU=response_london_bz.inp
       ../SCRIPTS/run.sh --outcmo --get="PAMXVC AOPROPER property_info"
       export MENU=vis_shield_natural_zz.inp
       ../SCRIPTS/run.sh --incmo --copy="PAMXVC AOPROPER property_info"
       rm -f PAMXVC AOPROPER property_info
   
       export MENU=response_london_symcon_bz.inp
       ../SCRIPTS/run.sh --outcmo --get="PAMXVC AOPROPER property_info"
       export MENU=vis_shield_symmetric_zz.inp
       ../SCRIPTS/run.sh --incmo  --copy="PAMXVC AOPROPER property_info"
       export MENU=vis_shield_none_zz.inp
       ../SCRIPTS/run.sh --incmo --copy="PAMXVC AOPROPER property_info"
   
   done
else
   for molecule in co2_*.mol; do
       export MOL=${molecule%\.mol}
   
       export MENU=response_london_bz.inp
       ../SCRIPTS/run.sh --outcmo --get="PAMXVC AOPROPER property_info"
       export MENU=vis_shield_natural_zz.inp
       ../SCRIPTS/run.sh --incmo  --copy="PAMXVC AOPROPER property_info"
       rm -f PAMXVC AOPROPER property_info
   
       export MENU=response_london_symcon_bz.inp
       ../SCRIPTS/run.sh --outcmo --get="PAMXVC AOPROPER property_info"
       export MENU=vis_shield_symmetric_zz.inp
       ../SCRIPTS/run.sh --incmo  --copy="PAMXVC AOPROPER property_info"
       export MENU=vis_shield_none_zz.inp
       ../SCRIPTS/run.sh --incmo --copy="PAMXVC AOPROPER property_info"
   
   done
fi

exit 0
