#!/bin/sh
#
# Set up calculation
# Since this calculation takes a long time to un we will only do this
# if the user asks for a benchmark
#
if [ ! -z "$RUNTYPE" ]; then
   case "${RUNTYPE}" in
      "small_benchmark"|"benchmark"|"large_benchmark")
         MENU=cc
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
MENU=cc;
export MENU;
###############################
echo "------>"
echo "runs using full symmetry"
echo
###############################
echo "C2v symmetry: eclipsed 1"
MOL=C2H4Cl2_ec1_c2v
export MOL;
../run.sh
###############################
echo "C2h symmetry: staggered"
MOL=C2H4Cl2_sta_c2h
export MOL;
../run.sh
###############################
echo "C2 symmetry: eclipsed 2"
MOL=C2H4Cl2_ec2_c2
export MOL;
../run.sh
###############################
if [ "$RUNTYPE" = "small_benchmark" ]; then
   echo
   echo "Completed small_benchmark"
   echo "<-----"
   exit
fi
###############################
echo
echo "------>"
echo "runs in lower symmetry"
echo
###############################
echo "Cs symmetry 2: eclipsed 1"
MOL=C2H4Cl2_ec1_cs2
export MOL;
../run.sh
###############################
echo "Ci symmetry: staggered"
MOL=C2H4Cl2_sta_ci
export MOL;
../run.sh
###############################
echo "No symmetry: eclipsed 2"
MOL=C2H4Cl2_ec2_c1
export MOL;
../run.sh
###############################
if [ "$RUNTYPE" = "benchmark" ]; then
   echo
   echo "Completed benchmark"
   echo "<-----"
   exit
fi
###############################
echo "Cs symmetry 1: eclipsed 1"
MOL=C2H4Cl2_ec1_cs1
export MOL;
../run.sh
###############################
echo "Cs symmetry 3: staggered"
MOL=C2H4Cl2_sta_cs
export MOL;
../run.sh
###############################
echo "C2 symmetry 1: eclipsed 1"
MOL=C2H4Cl2_ec1_c2
export MOL;
../run.sh
###############################
echo "C2 symmetry 2: staggered"
MOL=C2H4Cl2_sta_c2
export MOL;
../run.sh
###############################
echo "No symmetry: eclipsed 1"
MOL=C2H4Cl2_ec1_c1
export MOL;
../run.sh
###############################
echo "No symmetry: staggered"
MOL=C2H4Cl2_sta_c1
export MOL;
../run.sh
###############################
echo
echo "Completed large_benchmark"
echo "<-----"
