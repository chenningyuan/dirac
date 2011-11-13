#!/bin/sh
#
# Set up calculation
#
export MOL=LiH-c1

export MENU=gs
../run.sh -outcmo
cp DFCOEF DFCOEF.GS
export MENU=ion
../run.sh -incmo -outcmo
mv DFCOEF DFCOEF.ION
mv DFCOEF.GS DFCOEF
export MENU=stex
../run.sh -incmo -copy DFCOEF.ION
rm DFCOEF.ION DFCOEF
