#!/bin/sh
#
# Set up calculation
# Since this calculation takes a long time to un we will only do this
# if the user asks for a benchmark
#
if [ ! -z "$RUNTYPE" ]; then
   case "${RUNTYPE}" in
      "small_benchmark"|"benchmark"|"large_benchmark")
         MOL=C2H4Cl2_ec1_c2v
         ;;
      *)
      echo "we do not run the benchmark with runtype $RUNTYPE" 
         exit
         ;;
   esac
else
   echo "we do not run the benchmark with if the runtype is not specified (should be benchmark)" 
   exit
fi
export MOL;
###############################
echo "Moltra with scheme 4"
MENU=moltra4
export MENU;
../run.sh
###############################
echo "Moltra with scheme 6"
MENU=moltra6
export MENU;
../run.sh
###############################
echo "Completed benchmark"
