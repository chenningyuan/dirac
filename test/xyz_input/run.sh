#!/bin/bash

export MOL=CH3F.xyz
export MENU=xyz.inp
../SCRIPTS/run.sh

export MOL=CH3F.mol
export MENU=mol.inp
../SCRIPTS/run.sh

exit 0
