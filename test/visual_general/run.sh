#!/bin/sh

if [ "$RUNTYPE" = "short" ]; then
   for molecule in ne_d2h.mol; do
       export MOL=${molecule%\.mol}
   
       export MENU=response_ez_ez.inp
       ../SCRIPTS/run.sh --outcmo --get=PAMXVC
       export MENU=visual_zo.inp
       ../SCRIPTS/run.sh --incmo  --copy=PAMXVC --get=plot.3d.cube
       mv plot.3d.cube ${molecule%\.mol}_zo.3d
       export MENU=visual_ez_ez.inp
       ../SCRIPTS/run.sh --incmo  --copy=PAMXVC
   
       export MENU=response_bz_bz.inp
       ../SCRIPTS/run.sh --outcmo --get=PAMXVC
       export MENU=visual_bz_bz.inp
       ../SCRIPTS/run.sh --incmo  --copy=PAMXVC
   done
else
   for molecule in ne*mol; do
       export MOL=${molecule%\.mol}
   
       export MENU=response_ez_ez.inp
       ../SCRIPTS/run.sh --outcmo --get=PAMXVC
       export MENU=visual_zo.inp
       ../SCRIPTS/run.sh --incmo  --copy=PAMXVC --get=plot.3d.cube
       mv plot.3d.cube ${molecule%\.mol}_zo.3d
       export MENU=visual_ez_ez.inp
       ../SCRIPTS/run.sh --incmo  --copy=PAMXVC
   
       export MENU=response_bz_bz.inp
       ../SCRIPTS/run.sh --outcmo --get=PAMXVC
       export MENU=visual_bz_bz.inp
       ../SCRIPTS/run.sh --incmo  --copy=PAMXVC
   done
fi

export MOL=co.mol

export MENU=response_1z_bz.inp
../SCRIPTS/run.sh --outcmo --get=PAMXVC
export MENU=visual_1z_bz.inp
../SCRIPTS/run.sh --incmo  --copy=PAMXVC
export MENU=visual_1z_bz_only_para.inp
../SCRIPTS/run.sh --incmo  --copy=PAMXVC

export MENU=response_1z_2z.inp
../SCRIPTS/run.sh --outcmo --get=PAMXVC
export MENU=visual_1z_2z.inp
../SCRIPTS/run.sh --incmo  --copy=PAMXVC

# clean up
rm -f DFCOEF* PAMXVC*

exit 0
