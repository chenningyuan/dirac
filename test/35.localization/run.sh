#!/bin/ksh
# Calculate carbon atom
pam -mw 60 -nw 60 -utcmo -rmwrk -mpi 8 C C
cp DFCOEF DFC1AT
# Calculate hydrogen atom
/home/u1/trondsau/Dirac/pam -mw 60 -nw 60 -utcmo -mpi 8 -rmwrk H H
cp DFCOEF DFH1AT
ln -f -s DFH1AT DFH2AT
ln -f -s DFH1AT DFH3AT
ln -f -s DFH1AT DFH4AT
# Calculate methane
/home/u1/trondsau/Dirac/pam -mw 60 -nw 60 -utcmo -rmwrk -mpi 8 CH4 CH4
cp DFCOEF cf.CH4
# Pipek-Mezey localization using Mulliken population analysis
/home/u1/trondsau/Dirac/pam -mw 60 -nw 60 -incmo -utcmo -rmwrk -mpi 8 -rmwrk CH4 CH4locmul
cp DFCOEF cf.CH4locmul
# Pipek-Mezey localization using projection analysis
cp cf.CH4 DFCOEF
/home/u1/trondsau/Dirac/pam -nw 60 -mw 60 -incmo -copy DFC1AT -copy DFH1AT -copy DFH2AT -copy DFH3AT -copy DFH4AT -utcmo -rmwrk CH4 CH4locprj
cp DFCOEF cf.CH4locprj
exit 0
