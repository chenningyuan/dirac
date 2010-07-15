#!/bin/bash

export MOL=ne_c1.mol

export MENU=cme_sf.inp
../run.sh

export MENU=cme_dc.inp
../run.sh

exit 0
