#!/bin/sh

export MOL=H2O_nosym

export MENU=nmr_LinK_Jengine
../run.sh

export MENU=nmr_hermit
../run.sh

export MOL=H2O_df
export MENU=nmr_LinK_df-Jengine
../run.sh

exit 0
