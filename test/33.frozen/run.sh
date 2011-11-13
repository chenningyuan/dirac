#
#!/bin/sh
#
# Run calculation for oxygen atom
#
MOL=O;
MENU=O;
export MOL;
export MENU;
../run.sh -utcmo
mv DFCOEF cf.O
ln -s cf.O DFPROJ
#
# Run calculation for water
#
MOL=H2O;
MENU=H2O;
export MOL;
export MENU;
../run.sh -copy DFPROJ 
rm DFPROJ cf.O
