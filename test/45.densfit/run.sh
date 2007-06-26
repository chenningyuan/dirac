#!/bin/sh
#
# Set up calculation
#  
for mol in Ar; do
   MOL=$mol; export MOL;
   for method in direct; do
      for hamiltonian in NR 2C_SF 2C_SO 2C_GT; do
         for scheme in Jengine+df                  \
                       Jengine+LinK+df             \
                       Jengine+df+dfvxc_LDA        \
                       Jengine+df_LDA              \
                       Jengine+df+dfvxc_BLYP       \
                       Jengine+df_BLYP             \
                       Jengine+LinK+df+dfvxc_B3LYP \
                       Jengine+LinK+df_B3LYP       \
                       Jengine+df+dfvxc_B3LYP      \
                       Jengine+df_B3LYP            \
         ; do
            MENU=$hamiltonian"_"$method"_"$scheme; export MENU;
            if [ $hamiltonian == "2C_SF" ]; then
               cat "NR_"$method"_"$scheme".inp" | sed s/.NONREL/".IOTC\n.SPINFREE"/ |sed s/"#"// > $MENU".inp"
            elif [ $hamiltonian == "2C_SO" ]; then
               cat "NR_"$method"_"$scheme".inp" | sed s/.NONREL/".IOTC"/ |sed s/"#"// > $MENU".inp"
            elif [ $hamiltonian == "2C_GT" ]; then
               cat "NR_"$method"_"$scheme".inp" | sed s/.NONREL/".IOTC\n.GAUNT"/ |sed s/"#"// > $MENU".inp"
            fi
            ../run.sh
         done
      done
   done
done
