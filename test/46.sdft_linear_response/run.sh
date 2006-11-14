#!/bin/sh

MOL=he; export MOL;

MENU=blyp_collinear; export MENU;
../run.sh

MENU=blyp_noncollinear; export MENU;
../run.sh

MENU=lda_collinear; export MENU;
../run.sh

MENU=lda_noncollinear; export MENU;
../run.sh
