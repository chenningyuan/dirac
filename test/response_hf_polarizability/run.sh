#!/bin/sh

export MOL=H2O

export MENU=dc
../run.sh

export MENU=bss
../run.sh

export MENU=bss+mfsso
../run.sh

export MENU=bss_rkb+mfsso
../run.sh

exit 0
