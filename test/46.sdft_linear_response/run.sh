#!/bin/sh

MOL=hcl;   export MOL;

MENU=lda;  export MENU;
../run.sh

MENU=blyp; export MENU;
../run.sh
