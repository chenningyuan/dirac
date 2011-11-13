#
#!/bin/sh
#
# Run calculation for oxygen atom
#
echo "Four-component relativistic (DC) calculations"
MOL=O;
MENU=O;
export MOL;
export MENU;
../run.sh -outcmo
mv DFCOEF cf.O
ln -sf cf.O DFPROJ
#
# Run calculation for water
#
MOL=H2O;
MENU=H2O;
export MOL;
export MENU;
../run.sh -copy DFPROJ 
rm DFPROJ cf.O


###########################################
# run 2comp. IOTC calculations without AMFI
###########################################
#
# Run calculation for oxygen atom
#
echo "Two-component relativistic (IOTC) calculations"
MOL=O;
MENU=O.2c_iotc_noamfi;
export MOL;
export MENU;
../run.sh -utcmo
mv DFCOEF cf.O
ln -sf cf.O DFPROJ
#
# Run calculation for water
#
MOL=H2O;
MENU=H2O.2c_iotc_noamfi;
export MOL;
export MENU;
../run.sh -copy DFPROJ 
rm DFPROJ cf.O
