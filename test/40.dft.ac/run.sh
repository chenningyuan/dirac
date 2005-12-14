#!/bin/sh
#
# Set up calculation
#  
MOL=Ne; export MOL;
MENU=PBE0gracLB94; export MENU;
../run.sh
MENU=GLLBsaopLBalpha; export MENU;
../run.sh
