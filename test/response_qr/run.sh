#!/bin/sh
#
# Set up calculation
#
export MENU=qrpa1
export MOL=LiH
../SCRIPTS/run.sh
export MENU=tpa1
export MOL=He
../SCRIPTS/run.sh
export MENU=excprp1
export MOL=LiH
../SCRIPTS/run.sh
