#!/bin/bash


if [ $# -eq 0 ]
then
  echo "give the number of the test you want to check"
else

export INP_file=cc_uf6.inp
export MOL_file=UF6.xyz

if [ "$1" -eq 1 ] || [ "$1" -gt 8 ]
then
  FOLD=scheme1
  echo ""
  echo " compare results in ${FOLD}  "
  echo "======================================================"
  cd $FOLD
    fname="${INP_file%.*}_${MOL_file%.*}.out"
    grep -n -A 0 -a "Final CCSD" ${fname}
    echo "-----------------------------"
    grep -n -A 0 -a "Final CCSD" ../result/${fname}
    echo ""
    grep -n -A 1 -a "Final L1amp" ${fname}
    echo "-----------------------------"
    grep -n -A 1 -a "Final L1amp" ../result/${fname}
    echo ""
    grep -n -A 5 -a "Dipole length: X" ${fname}
    echo "-----------------------------"
    grep -n -A 5 -a "Dipole length: X" ../result/${fname}
  cd ..
  echo "======================================================"
fi

if [ "$1" -eq 2 ] || [ "$1" -gt 8 ]
then
  FOLD=scheme2
  echo ""
  echo " compare results in ${FOLD}  "
  echo "======================================================"
  cd $FOLD
    fname="${INP_file%.*}_${MOL_file%.*}.out"
    grep -n -A 0 -a "Final CCSD" ${fname}
    echo "-----------------------------"
    grep -n -A 0 -a "Final CCSD" ../result/${fname}
    echo ""
    grep -n -A 1 -a "Final L1amp" ${fname}
    echo "-----------------------------"
    grep -n -A 1 -a "Final L1amp" ../result/${fname}
    echo ""
    grep -n -A 5 -a "Dipole length: X" ${fname}
    echo "-----------------------------"
    grep -n -A 5 -a "Dipole length: X" ../result/${fname}
  cd ..
  echo "======================================================"
fi

if [ "$1" -eq 3 ] || [ "$1" -gt 8 ]
then
  FOLD=scheme3
  echo ""
  echo " compare results in ${FOLD}  "
  echo "======================================================"
  cd $FOLD
    fname="${INP_file%.*}_${MOL_file%.*}.out"
    grep -n -A 0 -a "Final CCSD" ${fname}
    echo "-----------------------------"
    grep -n -A 0 -a "Final CCSD" ../result/${fname}
    echo ""
    grep -n -A 1 -a "Final L1amp" ${fname}
    echo "-----------------------------"
    grep -n -A 1 -a "Final L1amp" ../result/${fname}
    echo ""
    grep -n -A 5 -a "Dipole length: X" ${fname}
    echo "-----------------------------"
    grep -n -A 5 -a "Dipole length: X" ../result/${fname}
  cd ..
  echo "======================================================"
fi

if [ "$1" -eq 4 ] || [ "$1" -gt 8 ]
then
  FOLD=scheme4
  echo ""
  echo " compare results in ${FOLD}  "
  echo "======================================================"
  cd $FOLD
    fname="${INP_file%.*}_${MOL_file%.*}.out"
    grep -n -A 0 -a "Final CCSD" ${fname}
    echo "-----------------------------"
    grep -n -A 0 -a "Final CCSD" ../result/${fname}
    echo ""
    grep -n -A 1 -a "Final L1amp" ${fname}
    echo "-----------------------------"
    grep -n -A 1 -a "Final L1amp" ../result/${fname}
    echo ""
    grep -n -A 5 -a "Dipole length: X" ${fname}
    echo "-----------------------------"
    grep -n -A 5 -a "Dipole length: X" ../result/${fname}
  cd ..
  echo "======================================================"
fi

if [ "$1" -eq 5 ] || [ "$1" -gt 8 ]
then
  FOLD=ccd
  echo ""
  echo " compare results in ${FOLD}  "
  echo "======================================================"
  cd $FOLD
    fname="${INP_file%.*}_${MOL_file%.*}.out"
    grep -n -A 0 -a "Final CCD" ${fname}
    echo "-----------------------------"
    grep -n -A 0 -a "Final CCD" ../result/${fname}
    echo ""
  cd ..
  echo "======================================================"
fi

export INP_FILE=cc_laf3.inp
export MOL_XYZ=LaF3.xyz

if [ "$1" -eq 6 ] || [ "$1" -gt 8 ]
then
  FOLD=script_exatensor
  echo ""
  echo " compare results in ${FOLD}  "
  echo "======================================================"
  cd $FOLD
    fname="${INP_file%.*}_${MOL_file%.*}.out"
    grep -n -A 0 -a "Final CCSD" ${fname}
    echo "-----------------------------"
    grep -n -A 0 -a "Final CCSD" ../result/${fname}
    echo ""
    grep -n -A 1 -a "Final L1amp" ${fname}
    echo "-----------------------------"
    grep -n -A 1 -a "Final L1amp" ../result/${fname}
    echo ""
    grep -n -A 5 -a "Dipole length: X" ${fname}
    echo "-----------------------------"
    grep -n -A 5 -a "Dipole length: X" ../result/${fname}
  cd ..
  echo "======================================================"
fi

export INP_file=cc_sf6.inp
export MOL_file=SF6.xyz
if [ "$1" -eq 7 ] || [ "$1" -gt 8 ]
then
  FOLD=script_talsh
  echo ""
  echo " compare results in ${FOLD}  "
  echo "======================================================"
  cd $FOLD
    fname="${INP_file%.*}_${MOL_file%.*}.out"
    grep -n -A 0 -a "Final CCSD" ${fname}
    echo "-----------------------------"
    grep -n -A 0 -a "Final CCSD" ../result/${fname}
    echo ""
  cd ..
  echo "======================================================"
fi

fi
