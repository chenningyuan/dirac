#!/bin/sh

MOL=H2O2;          export MOL;

MENU=shielding_dc; export MENU;
../run.sh

MENU=shielding_nr; export MENU;
../run.sh

MENU=spinspin_dc;  export MENU;
../run.sh

MENU=spinspin_nr;  export MENU;
../run.sh

exit 0
