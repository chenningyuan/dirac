#!/bin/sh
#
# Set up calculation
#  
MENU=PBE0gracLB94; export MENU;
MOL=Ne; export MOL;
../run.sh
MENU=GLLBsaopLBalpha; export MENU;
../run.sh
MOL=He; export MOL;
../run.sh
MOL=Be; export MOL;
../run.sh
MOL=Ar; export MOL;
../run.sh
