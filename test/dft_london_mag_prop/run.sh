#!/bin/sh

export MOL=hehe_c1
export MENU=lda_lvcorr_sdft_london_shield.inp
../run.sh
export MENU=lda_lvcorr_sdft_london_shield_gauge007.inp
../run.sh
export MENU=lda_lvcorr_nosdft_london_shield.inp
../run.sh
export MENU=lda_lvcorr_nosdft_london_shield_gauge007.inp
../run.sh


#
#
#export MENU=lda_lvcorr_sdft_london_shield.inp
#../run.sh -put grid_file
#export MENU=lda_lvcorr_sdft_london_shield_gauge007.inp
#../run.sh -put grid_file
#export MENU=lda_lvcorr_nosdft_london_shield.inp
#../run.sh -put grid_file
#export MENU=lda_lvcorr_nosdft_london_shield_gauge007.inp
#../run.sh -put grid_file

##OLD
#
#export MOL=hehe_c1
#export MENU=lda_levyle_nosdft_london_shield.inp
#../run.sh
##export MENU=blyp_levyle_nosdft_london_shield.inp
##../run.sh
##export MENU=b3lyp_levyle_nosdft_london_shield.inp
##../run.sh
#export MENU=lda_lvcorr_nosdft_london_shield.inp
#../run.sh
#
#
#export MOL=hehe_c1
#export MENU=lda_levyle_nosdft_london_shield_gauge007.inp
#../run.sh
#export MENU=lda_lvcorr_nosdft_london_shield_gauge007.inp
#../run.sh
##export MOL=hehe_c1_shifted
##export MENU=lda_levyle_nosdft_london_shield.inp
##../run.sh
##export MENU=blyp_levyle_nosdft_london_shield.inp
##../run.sh
##export MENU=b3lyp_levyle_nosdft_london_shield.inp
##../run.sh
##export MENU=lda_lvcorr_nosdft_london_shield.inp
##../run.sh
#
### SDFT
#export MOL=hehe_c1
#export MENU=lda_levyle_sdft_london_shield.inp
#../run.sh
#export MENU=blyp_levyle_sdft_london_shield.inp
##../run.sh
#export MENU=b3lyp_levyle_sdft_london_shield.inp
##../run.sh
#export MENU=lda_levyle_sdft_london_shield_gauge007.inp
#../run.sh
#export MENU=lda_lvcorr_sdft_london_shield.inp
#../run.sh
#export MENU=lda_lvcorr_sdft_london_shield_gauge007.inp
#../run.sh
#
#export MOL=hehe_c1_shifted
#export MENU=lda_levyle_sdft_london_shield.inp
##../run.sh
#export MENU=blyp_levyle_sdft_london_shield.inp
##../run.sh
#export MENU=b3lyp_levyle_sdft_london_shield.inp
##../run.sh
#export MENU=lda_lvcorr_sdft_london_shield.inp
##../run.sh
#
#
#
## SDFT
## now tests are compared with results from table II and table IV from 
## J. Chem. Phys. 128, 104101 (2008)
## there 'uncoupled' and 'coupled' mDKS-RMB
## can be compared with our 'dft-nosdft' and 'dft-sdft' respectively
#
#export MOL=hf
#export MENU=lda_levyle_nosdft_london_shield.inp
##../run.sh -copy grid_file
#
#export MOL=hf_shifted
#export MENU=lda_levyle_nosdft_london_shield.inp
##../run.sh
#
#
#
#export MOL=hehe_c1
#export MENU=lda_levyle_sdft_london_shield.inp
##../run.sh
##../run.sh -copy grid_file
#
#
#export MOL=hf_shifted
#export MENU=lda_levyle_sdft_london_shield.inp
##../run.sh
#
#
#
#
exit 0
