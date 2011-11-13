#!/bin/sh
#
# Run calculation for Hydrogen
#
MOL=H;
MENU=H;
export MOL;
export MENU;
../run.sh -utcmo
mv DFCOEF cf.H
ln -s cf.H DFHAT1
ln -s cf.H DFHAT2
ln -s cf.H DFHAT3
ln -s cf.H DFHAT4
#
# Run calculation for Carbon
#
MOL=C;
MENU=C;
export MOL;
export MENU;
../run.sh -utcmo
mv DFCOEF cf.C
ln -s cf.C DFCAT
#
# Run calculation for Methane
#
MOL=CH4;
MENU=CH4;
export MOL;
export MENU;
../run.sh -copy DFCAT -copy DFHAT1 -copy DFHAT2 -copy DFHAT3 -copy DFHAT4
rm DFCAT DFHAT1 DFHAT2 DFHAT3 DFHAT4 cf.H cf.C
