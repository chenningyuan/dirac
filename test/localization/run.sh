#
#!/bin/sh
#
# Run calculation for carbon atom
#
MOL=C;
MENU=C;
export MOL;
export MENU;
../run.sh -outcmo
mv DFCOEF cf.C
ln -f -s cf.C DFC1AT
#
# Run calculation for hydrogen
#
MOL=H;
MENU=H;
export MOL;
export MENU;
../run.sh -outcmo
mv DFCOEF cf.H
ln -f -s cf.H DFH1AT
ln -f -s cf.H DFH2AT
ln -f -s cf.H DFH3AT
ln -f -s cf.H DFH4AT
#
# Run calculation for methane
#
MOL=CH4;
MENU=CH4;
export MOL;
export MENU;
../run.sh -outcmo
#
# Pipek-Mezey localization using Mulliken population analysis
#
MOL=CH4;
MENU=CH4locmul;
export MOL;
export MENU;
../run.sh -incmo
#
# Pipek-Mezey localization using projection analysis
#
MOL=CH4;
MENU=CH4locprj;
export MOL;
export MENU;
../run.sh -incmo -copy DFC1AT -copy DFH1AT -copy DFH2AT -copy DFH3AT -copy DFH4AT
rm cf.C cf.H DFCOEF DFC1AT DFH1AT DFH2AT DFH3AT DFH4AT
exit 0
