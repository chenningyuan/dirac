#/bin/sh

# This script is based on the implementation in Ref. https://arxiv.org/pdf/1701.07862
# It runs the program AVAS, that computes overlap between a big basis set B1 of a molecule
# and a minimal basis B2 (containing for instance the 3d orbitals of the metallic center of the molecule).
# In this script, the Intrinsic Atomic Orbitals are first constructed. The valence IAO of the fragment (B2frag)
# are then considered to form B2.

#################
# MANUAL SET-UP #
#################
path_to_geniao=../../build-caar-ibos/ # PATH TO THE EXECUTABLE
path_to_avas=../../build-caar-ibos/ # PATH TO THE EXECUTABLE
inp1=HF_closedshell # INPUT FOR THE MOLECULE (ionic species can be considered)
inp2=HF # INPUT FOR THE FRAGMENT (neutral)
mol=HeNe # NAME OF THE MOLECULE
charge=2 # CHARGE OF THE MOLECULE 
nelec=10 # TOTAL NUMBER OF ELECTRONS
basis1_iao=sto-3g # BASIS SET OF THE ORIGINAL PROBLEM
basis2_iao=cc-pvtz # LARGE BASIS SET FOR THE FRAGMENTS
nfrag=2 # NUMBER OF FRAGMENTS (atoms for now!)
# FRAGMENTS NAMES, ATOMIC NUMBER, COORDINATES
frag[1]=He
Z[1]=002 # MUST BE IN THIS FORMAT
coord[1]="0. 0. -1.0"
frag[2]=Ne
Z[2]=010 # MUST BE IN THIS FORMAT
coord[2]="0. 0. +1.0"

# Choose the atomic valence orbitals B2.
B2frag=Ne # atom for basis B2
atom_number=010 # The atomic number of B2frag (format i3.3) that will form the basis B2 (only the valence orbitals will then be considered by the program)
basis_avas=cc-pvtz

######################
# GENERATE XYZ FILES #
######################
cat > MOLECULE.XYZ << EOF
$nfrag
molecule
EOF
for i in `seq 1 $nfrag` ; do
echo ${frag[$i]} ${coord[$i]} >> MOLECULE.XYZ
cat > ${frag[$i]}.xyz << EOF
1
fragment
${frag[$i]} 0. 0. 0.
EOF
done
cp MOLECULE.XYZ $mol.xyz

############################################
# GENERATE INPUT FILE (charge and version) #
############################################
echo "Simple_2013" > INPUT # choices: Simple_2014, Simple_2013, Stndrd_2013
echo "$charge" >> INPUT

echo "$nelec" > INPUT_AVAS
echo "$charge" >> INPUT_AVAS
echo "$atom_number" >> INPUT_AVAS

cat > $inp1.inp << EOF
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

cat > $inp2.inp << EOF
**DIRAC
.WAVE FUNCTION
**WAVE FUNCTION
.SCF
*SCF
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
# for IAO
pam --replace basis=$basis1_iao --replace nelec=$nelec --inp=$inp1.inp --mol=$mol.xyz --get="DFCOEF" --silent
mv DFCOEF MOLECULE.dfcoef
mv "$inp1"_"$mol"_basis="$basis1_iao"_nelec="$nelec".out "$mol"_"$basis1_iao"_nelec"$nelec".out

# for AVAS
pam --replace basis=$basis_avas --replace nelec=$nelec --inp=$inp1.inp --mol=$mol.xyz --get="DFCOEF" --silent
mv DFCOEF B1.dfcoef
mv "$inp1"_"$mol"_basis="$basis_avas"_nelec="$nelec".out "$mol"_"$basis_avas"_nelec"$nelec".out

for i in `seq 1 $nfrag` ; do
  pam --replace basis=$basis2_iao --inp=$inp2.inp --mol=${frag[$i]}.xyz --get="DFCOEF" --silent
  mv DFCOEF ${Z[$i]}.dfcoef
  mv "$inp2"_"${frag[$i]}"_basis="$basis2_iao".out "$mol"_"$basis2_iao".out
done

echo
echo
echo "#################"
echo "# GENERATE IAOs #"
echo "#################"
echo

$path_to_geniao/geniao.x

echo
echo
echo "############################################################################"
echo "# TEST IF MY TRANSFORMATION TO GENERATE THE INITIAL DFCOEF FILE IS CORRECT #"
echo "############################################################################"
echo
cp test.dfcoef DFCOEF
pam --replace basis=$basis1_iao --replace nelec=$nelec --inp=$inp1.inp --mol=$mol.xyz --put="DFCOEF" --silent
mv "$inp1"_"$mol"_basis="$basis1_iao"_nelec="$nelec".out "$mol"_"$basis1_iao"_nelec"$nelec"_restarted.out
convergence=10
convergence=$(grep "* Convergence after" "$mol"_"$basis1_iao"_nelec"$nelec"_restarted.out | awk '{print $4}')
if (( $(echo "$convergence != 1" | bc -l) )); then
  echo
  echo
  echo "Convergence after $convergence iterations"
  echo ERROR IN THE RECONSTRUCTION OF THE DFCOEF FILE
else
  echo
  echo
  echo Convergence after $convergence iterations
  echo RECONSTRUCTION OF THE DFCOEF FILE SUCCEEDED
fi

echo
echo
echo "########################################################"
echo "# RESTART HF CALCULATION WITH THE NEW IAO COEFFICIENTS #"
echo "########################################################"
echo
cp iao.dfcoef DFCOEF
pam --replace basis=$basis1_iao --replace nelec=$nelec --inp=$inp1.inp --mol=$mol.xyz --put="DFCOEF" --silent
mv "$inp1"_"$mol"_basis="$basis1_iao"_nelec="$nelec".out "$mol"_"$basis1_iao"_nelec"$nelec"_iao.out

convergence=10
convergence=$(grep "* Convergence after" "$mol"_"$basis1_iao"_nelec"$nelec"_iao.out | awk '{print $4}')
iteration1=$(grep "It.    1" "$mol"_"$basis1_iao"_nelec"$nelec".out | tail -1 | awk '{print $3}')
iteration1_iao=$(grep "It.    1" "$mol"_"$basis1_iao"_nelec"$nelec"_iao.out | tail -1 | awk '{print $3}')
HF_energy=$(grep "   Total energy                             :" "$mol"_"$basis1_iao"_nelec"$nelec".out | awk '{print $4}')

if (( $(echo "$convergence != 1" | bc -l) )); then
  echo
  echo
  echo "Convergence after $convergence iterations"
  echo "Energy of the first calculation at iteration 1  :" $iteration1
  echo "Energy at iteration 1 using the IAO coeffcients :" $iteration1_iao
  echo "Hartree-Fock energy                             :" $HF_energy
  echo ERROR IN THE GENERATION OF IAO COEFFICIENTS
else
  echo
  echo
  echo Convergence after $convergence iterations
  echo Hartree-Fock energy: $HF_energy
  echo GENERATION OF IAO COEFFICIENTS CORRECTLY SUCCEEDED
fi

echo 
echo 
echo "#############################"
echo "# EXECUTE THE AVAS PROGRAM  #"
echo "#############################"
echo 
cp "$atom_number"_iao.dfcoef atom.dfcoef
$path_to_avas/avas.x

echo
echo
echo "################################################################"
echo "# RESTART HF CALCULATION WITH THE NEW TRANSFORMED COEFFICIENTS #"
echo "################################################################"
echo
cp avas.dfcoef DFCOEF
pam --replace basis=$basis_avas --replace nelec=$nelec --inp=$inp1.inp --mol=$mol.xyz --put="DFCOEF" --silent
mv "$inp1"_"$mol"_basis="$basis_avas"_nelec="$nelec".out "$mol"_"$basis_avas"_nelec"$nelec"_avas.out

convergence=10
convergence=$(grep "* Convergence after" "$mol"_"$basis_avas"_nelec"$nelec"_avas.out | awk '{print $4}')
iteration1=$(grep "It.    1" "$mol"_"$basis_avas"_nelec"$nelec".out | tail -1 | awk '{print $3}')
iteration1_avas=$(grep "It.    1" "$mol"_"$basis_avas"_nelec"$nelec"_avas.out | tail -1 | awk '{print $3}')
HF_energy=$(grep "   Total energy                             :" "$mol"_"$basis_avas"_nelec"$nelec".out | awk '{print $4}')

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
