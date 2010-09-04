#!/bin/sh
#
# Run calculation for Hydrogen in linear symmetry
#
MOL=H_lin;
MENU=H_lin;
export MOL;
export MENU;
../run.sh -get DFACMO
#
# Export coefficients to Cs symmetry
#
MOL=H_Cs;
MENU=H_Cs;
export MOL;
export MENU;
../run.sh -copy DFACMO -outcmo 
mv DFCOEF cf.H
rm DFACMO
ln -sf cf.H DFHAT1
ln -sf cf.H DFHAT2
ln -sf cf.H DFHAT3
#
# Run calculation for Nitrogen in linear symmetry
#
MOL=N_lin;
#MENU=N_lin;
MENU=N_lin.noamfi;
export MOL;
export MENU;
../run.sh -get DFACMO
mv DFCOEF cf.C
ln -s cf.C DFCAT
#
# Export coefficients to Cs symmetry
#
MOL=N_Cs;
#MENU=N_Cs;
MENU=N_Cs.noamfi;
export MOL;
export MENU;
../run.sh -copy DFACMO -outcmo 
mv DFCOEF cf.N
rm DFACMO
ln -sf cf.N DFNNNN
#
# Run calculation for Ammonia
#
MOL=NH3;
#MENU=NH3;
MENU=NH3.noamfi;
export MOL;
export MENU;
../run.sh -copy DFNNNN -copy DFHAT1 -copy DFHAT2 -copy DFHAT3
rm DFNNNN DFHAT1 DFHAT2 DFHAT3 cf.H cf.N
