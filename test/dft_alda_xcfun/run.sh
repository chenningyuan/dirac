#!/bin/sh

export MOL=he

export MENU=b3lyp.inp
../run.sh
export MENU=b3lyp_alda.inp
../run.sh
export MENU=b3lyp_xalda.inp
../run.sh

export MOL=zn

export MENU=pbe.inp
../run.sh
export MENU=pbe_alda+.inp
../run.sh
export MENU=pbe_alda-.inp
../run.sh

exit 0
