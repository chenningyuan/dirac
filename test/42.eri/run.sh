#!/bin/sh
#
# Set up calculation
#  
for code in eri hermit; do
   MOL="H2O2_"$code; export MOL;
   for method in direct; do
      MENU="NR_"$method"_"$code; export MENU;
      ../run.sh
   done
done
