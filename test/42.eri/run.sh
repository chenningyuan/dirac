#!/bin/sh
#
# Set up calculation
#  
for mol in H2O2_eri; do
   MOL=$mol; export MOL;
   for method in direct; do
     # for hamiltonian in NR 2C_SF 2C_noamfi 2C_SO 2C_GT; do
      for hamiltonian in NR 2C_SF 2C_noamfi; do
         MENU=$hamiltonian"_"$method"_eri"; export MENU;
         ../run.sh
      done
   done
done
