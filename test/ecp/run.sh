#!/bin/bash

export MENU=ecp_d2h

export MOL=f_arpp_d2h
../SCRIPTS/run.sh

export MOL=f_sopp_d2h
../SCRIPTS/run.sh

exit 0
