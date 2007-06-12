#!/bin/sh

MOL=he;              export MOL;

MENU=blyp;           export MENU;
../run.sh
MENU=blyp_collinear; export MENU;
../run.sh
MENU=blyp_alda;      export MENU;
../run.sh

MENU=lda;            export MENU;
../run.sh
MENU=lda_collinear;  export MENU;
../run.sh
