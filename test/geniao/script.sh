#/bin/sh

#################
# MANUAL SET-UP #
#################
path_to_geniao=../../build-caar-ibos/ # PATH TO THE EXECUTABLE
inp1=HF_closedshell # INPUT FOR THE MOLECULE (ionic species can be considered)
inp2=HF # INPUT FOR THE FRAGMENT (neutral)
mol=HeNe # NAME OF THE MOLECULE
charge=2 # CHARGE OF THE MOLECULE 
nelec=10 # TOTAL NUMBER OF ELECTRONS
basis1=sto-3g # BASIS SET OF THE ORIGINAL PROBLEM
basis2=cc-pvtz # LARGE BASIS SET FOR THE FRAGMENTS
nfrag=2 # NUMBER OF FRAGMENTS (atoms for now!)
# FRAGMENTS NAMES, ATOMIC NUMBER, COORDINATES
frag[1]=He 
Z[1]=002 # MUST BE IN THIS FORMAT
coord[1]="0. 0. -1."
frag[2]=Ne
Z[2]=010
coord[2]="0. 0. +1."

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
echo "Simple_2014" > INPUT # choices: Simple_2014, Simple_2013, Stndrd_2013
echo "$charge" >> INPUT

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

#############################
# GENERATE THE DFCOEF FILES #
#############################
pam --replace basis=$basis1 --replace nelec=$nelec --inp=$inp1.inp --mol=$mol.xyz --get="DFCOEF" --silent
mv DFCOEF MOLECULE.dfcoef
mv "$inp1"_"$mol"_basis="$basis1"_nelec="$nelec".out "$mol"_"$basis1"_nelec"$nelec".out

for i in `seq 1 $nfrag` ; do
  pam --replace basis=$basis2 --inp=$inp2.inp --mol=${frag[$i]}.xyz --get="DFCOEF" --silent
  mv DFCOEF ${Z[$i]}.dfcoef
  mv "$inp2"_"$mol"_basis="$basis2".out "$mol"_"$basis2".out
done

#################
# GENERATE IAOs #
#################
$path_to_geniao/geniao.x

########################################################
# RESTART HF CALCULATION WITH THE NEW IAO COEFFICIENTS #
########################################################
cp iao.dfcoef DFCOEF
pam --replace basis=$basis1 --replace nelec=$nelec --inp=$inp1.inp --mol=$mol.xyz --put="DFCOEF" --silent
mv "$inp1"_"$mol"_basis="$basis1"_nelec="$nelec".out "$mol"_"$basis1"_nelec"$nelec"_restarted.out

convergence=10
convergence=$(grep "* Convergence after" "$mol"_"$basis1"_nelec"$nelec"_restarted.out | awk '{print $4}')

if (( $(echo "$convergence != 1" | bc -l) )); then
  echo
  echo
  echo Convergence after $convergence iterations
  echo ERROR IN THE GENERATION OF IAO COEFFICIENTS
else
  echo
  echo
  echo Convergence after $convergence iterations
  echo GENERATION OF IAO COEFFICIENTS CORRECTLY SUCCEEDED
fi
