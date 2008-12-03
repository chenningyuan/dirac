#!/bin/sh

export MOL=hf

export MENU=b3lyp.inp
../run.sh -copy d-aug-cc-pVDZ

export MENU=blyp.inp
../run.sh -copy d-aug-cc-pVDZ

export MENU=lda.inp
../run.sh -copy d-aug-cc-pVDZ

exit 0
