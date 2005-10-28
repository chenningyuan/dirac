#!/bin/sh
#
# Set up calculation
#
export MENU=qrpa1
export MOL=LiH
../run.sh
export MENU=tpa1
export MOL=He
../run.sh
export MENU=excprp1
export MOL=LiH
../run.sh
