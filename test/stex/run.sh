#!/bin/sh
#
# Set up calculation
#
export MOL=LiH-c1

export MENU=gs
../SCRIPTS/run.sh -outcmo
cp DFCOEF DFCOEF.GS
export MENU=ion
../SCRIPTS/run.sh -incmo -outcmo
mv DFCOEF DFCOEF.ION
mv DFCOEF.GS DFCOEF
export MENU=stex
../SCRIPTS/run.sh -incmo -copy DFCOEF.ION
rm DFCOEF.ION DFCOEF
