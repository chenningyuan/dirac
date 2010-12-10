#!/bin/sh

export MOL=H2O

export MENU=dc
../SCRIPTS/run.sh

export MENU=bss
../SCRIPTS/run.sh

export MENU=bss+mfsso
../SCRIPTS/run.sh

export MENU=bss_rkb+mfsso
../SCRIPTS/run.sh

exit 0
