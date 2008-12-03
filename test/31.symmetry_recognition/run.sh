#!/bin/sh

export MENU=lda

export MOL=H2O
../run.sh

export MOL=CH4
../run.sh

export MOL=SF6
../run.sh

exit 0
