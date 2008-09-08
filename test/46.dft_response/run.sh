#!/bin/sh

MOL=hf;         export MOL;

MENU=b3lyp.inp; export MENU;
../run.sh -copy d-aug-cc-pVDZ

MENU=blyp.inp;  export MENU;
../run.sh -copy d-aug-cc-pVDZ

MENU=lda.inp;   export MENU;
../run.sh -copy d-aug-cc-pVDZ

exit 0
