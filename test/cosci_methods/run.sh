#!/bin/sh
#
# Set up calculation
#

MOL=O;
export MOL;
###############################
echo "Levy-Leblond non-relativistic Hamiltonian "
MENU=ci.1;
export MENU;
../run.sh
###############################
echo "Dirac-Coulomb spin-free Hamiltonian "
MENU=ci.2;
export MENU;
../run.sh -utcmo
###############################
###############################
if [ "$RUNTYPE" = "short" ]; then
   echo
   echo "Completed short test"
   echo "<-----"
   exit
fi
###############################
echo "Dirac-Coulomb Hamiltonian using spin-free MO in CI"
MENU=ci.3;
export MENU;
../run.sh -incmo
###############################
echo "Full Dirac-Coulomb Hamiltonian"
MENU=ci.4;
export MENU;
../run.sh
rm DFCOEF
################################
### TWO COMPONENT CALCULATIONS #
################################
echo "Barysz-Sadlej-Snijders spin-free (from beginning, BSSsfB) Hamiltonian "
MENU=O.bss_sfB.scf_ci;
export MENU;
../run.sh -utcmo
###############################
#REMOVE_BEGIN
echo "BSS+MFSSO Hamiltonian at post-SCF"
MENU=O.bss+mfsso.postscf.ci;
export MENU;
../run.sh -incmo
rm DFCOEF
###############################
echo "Two-component BSS+MFSSO Hamiltonian "
MENU=O.bss+mfsso.scf_ci;
export MENU;
../run.sh 
#REMOVE_END
###############################
echo "Two-component BSS Hamiltonian "
MENU=O.bss.scf_ci;
export MENU;
../run.sh 
###############################

