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
if [ "$RUNTYPE" = "short" ]; then
   echo
   echo "Completed short test"
   echo "<-----"
   exit
fi
MOL=Be; export MOL;
../run.sh
MOL=Ar; export MOL;
../run.sh
