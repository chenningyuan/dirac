#/bin/sh

path_to_geniao=../../../build-caar-ibos/
inp=HF_closedshell
mol=He2
charge=0
nelec=4
basis1=sto-3g
basis2=cc-pvtz

cp $mol.xyz MOLECULE.XYZ
echo "Simple_2014" > INPUT # choices: Simple_2014, Simple_2013, Stndrd_2013
echo "$charge" >> INPUT

pam --replace basis=$basis1 --replace nelec=$nelec --inp=$inp.inp --mol=$mol.xyz --get="DFCOEF" --silent
mv DFCOEF MOLECULE.dfcoef
mv "$inp"_"$mol"_basis="$basis1"_nelec="$nelec".out "$mol"_"$basis1"_nelec"$nelec".out

inp=HF
mol=He
pam --replace basis=$basis2 --inp=$inp.inp --mol=$mol.xyz --get="DFCOEF" --silent
mv DFCOEF 002.dfcoef
mv "$inp"_"$mol"_basis="$basis2".out "$mol"_"$basis2".out

$path_to_geniao/geniao.x

cp iao.dfcoef DFCOEF
inp=HF_closedshell
mol=He2
pam --replace basis=$basis1 --replace nelec=$nelec --inp=$inp.inp --mol=$mol.xyz --put="DFCOEF" --silent
mv "$inp"_"$mol"_basis="$basis1"_nelec="$nelec".out "$mol"_"$basis1"_nelec"$nelec"_restarted.out

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

/bin/rm *dfcoef *dfpcmo DFCOEF
