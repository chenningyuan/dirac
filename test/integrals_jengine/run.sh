#!/bin/sh
#
# Set up calculation
#  
for mol in H2O2; do
   MOL=$mol; export MOL;
   for method in direct; do
      for hamiltonian in NR 2C_SF 2C_SO 2C_GT; do
         MENU=$hamiltonian"_"$method"_Jengine+LinK"; export MENU;
         ../SCRIPTS/run.sh
         MENU=$hamiltonian"_"$method"_Jengine+LinK_LDA"; export MENU;
         ../SCRIPTS/run.sh
      done
   done
done
