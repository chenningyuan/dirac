#/bin/sh

# This script is based on the implementation in Ref. https://arxiv.org/pdf/1701.07862
# It runs the program AVAS, that computes overlap between a big basis set B1 of a molecule
# and a minimal basis B2 (containing for instance the 3d orbitals of the metallic center of the molecule).
# In the program, the valence orbitals of the fragment (called B2frag) will be considered.

#################
# MANUAL SET-UP #
#################
path_to_avas=../../build-caar-ibos/ # PATH TO THE EXECUTABLE
inp=HF # INPUT FOR THE MOLECULE (ionic species can be considered)

#molecule:
mol=HeNe # NAME OF THE MOLECULE
nelec=10
charge=2
basis1=cc-pvtz # LARGE BASIS SET B1
cat > $mol.xyz << EOF
2
$mol
He 0. 0. -0.5
Ne 0. 0. +0.5
EOF

#atom:
B2frag=Ne # atom for basis B2
atom_number=010 # The atomic number of B2frag (format i3.3) that will form the basis B2 (only the valence orbitals will then be considered by the program)
nelec2=10
basis2=cc-pvdz # can be big, anyway only the valence orbitals are considered by the program !
cat > $B2frag.xyz << EOF
1
atom
$B2frag 0. 0. 0.
EOF

#######################
# GENERATE INPUT FILE #
#######################
echo "$nelec" > INPUT_AVAS
echo "$charge" >> INPUT_AVAS
echo "$atom_number" >> INPUT_AVAS

cat > $inp.inp << EOF
**DIRAC
.WAVE FUNCTION
**WAVE FUNCTION
.SCF
*SCF
.CLOSED SHELL
nelec
.MAXITR
1000
**HAMILTONIAN
.NONREL
**GENERAL
.PCMOUT
**INTEGRALS
.NUCMOD
 1
*READIN
.UNCONTRACT
**MOLECULE
*SYMMETRY
.NOSYM
*BASIS
.DEFAULT 
basis
*END OF INPUT
EOF

echo 
echo 
echo "#############################"
echo "# GENERATE THE DFCOEF FILES #"
echo "#############################"
pam --replace basis=$basis1 --replace nelec=$nelec --inp=$inp.inp --mol=$mol.xyz --get="DFCOEF" --silent
mv DFCOEF B1.dfcoef
mv "$inp"_"$mol"_basis="$basis1"_nelec="$nelec".out "$mol"_"$basis1"_nelec"$nelec".out

pam --replace basis=$basis2 --replace nelec=$nelec2 --inp=$inp.inp --mol=$B2frag.xyz --get="DFCOEF" --silent
mv DFCOEF atom.dfcoef
mv "$inp"_"$B2frag"_basis="$basis2"_nelec="$nelec2".out "$B2frag"_"$basis2"_nelec"$nelec2".out

echo 
echo 
echo "#############################"
echo "# EXECUTE THE AVAS PROGRAM  #"
echo "#############################"
$path_to_avas/avas.x

echo
echo
echo "################################################################"
echo "# RESTART HF CALCULATION WITH THE NEW TRANSFORMED COEFFICIENTS #"
echo "################################################################"
echo
cp avas.dfcoef DFCOEF
pam --replace basis=$basis1 --replace nelec=$nelec --inp=$inp.inp --mol=$mol.xyz --put="DFCOEF" --silent
mv "$inp"_"$mol"_basis="$basis1"_nelec="$nelec".out "$mol"_"$basis1"_nelec"$nelec"_avas.out

convergence=10
convergence=$(grep "* Convergence after" "$mol"_"$basis1"_nelec"$nelec"_avas.out | awk '{print $4}')
iteration1=$(grep "It.    1" "$mol"_"$basis1"_nelec"$nelec".out | tail -1 | awk '{print $3}')
iteration1_avas=$(grep "It.    1" "$mol"_"$basis1"_nelec"$nelec"_avas.out | tail -1 | awk '{print $3}')
HF_energy=$(grep "   Total energy                             :" "$mol"_"$basis1"_nelec"$nelec".out | awk '{print $4}')

if (( $(echo "$convergence != 1" | bc -l) )); then
  echo
  echo
  echo "Convergence after $convergence iterations"
  echo "Energy of the first calculation at iteration 1     :" $iteration1
  echo "Energy at iteration 1 with transformed coeffcients :" $iteration1_avas
  echo "Hartree-Fock energy                                :" $HF_energy
  echo ERROR IN THE TRANSFORMATION, NOT UNITARY ?
else
  echo
  echo
  echo Convergence after $convergence iterations
  echo Hartree-Fock energy: $HF_energy
  echo TRANSFORMED ORBITALS GIVE THE CORRECT ENERGY
fi
