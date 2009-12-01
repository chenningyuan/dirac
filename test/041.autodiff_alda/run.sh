#!/bin/sh

export MOL=he

export MENU=b3lyp.inp
../run.sh
export MENU=b3lyp_alda.inp
../run.sh
export MENU=b3lyp_xalda.inp
../run.sh
export MENU=blyp_alda-.inp
../run.sh

exit 0
