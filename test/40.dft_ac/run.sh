#!/bin/sh

MOL=Ne;               export MOL;

MENU=PBE0gracLB94;    export MENU;
../run.sh

MENU=GLLBsaopLBalpha; export MENU;
../run.sh

exit 0
