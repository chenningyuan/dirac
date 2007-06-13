#!/bin/sh
#
# Set up calculation
#  
for code in eri Jengine+LinK hermit; do
   MOL="H2O2_"$code; export MOL;
   for method in direct; do
      for hamiltonian in NR 2C_SF 2C_SO 2C_GT; do
         MENU=$hamiltonian"_"$method"_"$code; export MENU;
         ../run.sh
      done
   done
done
