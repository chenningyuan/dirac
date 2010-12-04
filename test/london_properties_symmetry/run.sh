#!/bin/sh

# gosia: this test will be merged to dft_lao tests
# i will modify it, now it takes too long

export MOL=co2_c1.mol
export MENU=lda_levyle_sdft_shielding.inp
#../run.sh
export MENU=lda_levyle_nosdft_shielding.inp
#../run.sh
export MENU=lda_lvcorr_sdft_shielding.inp
#../run.sh
export MENU=lda_lvcorr_nosdft_shielding.inp
#../run.sh

#export MENU=lda_levyle_nosdft_shielding.inp
#../run.sh -get test_ar -get dft_before -get dft_after
#mv test_ar    test_ar.c1
#mv dft_before dft_before.c1
#mv dft_after  dft_after.c1
#
#
#export MOL=co2_d2h.mol
#export MENU=lda_levyle_nosdft_shielding.inp
#../run.sh -get test_ar -get dft_before -get dft_after
#mv test_ar test_ar.d2h
#mv dft_before dft_before.d2h
#mv dft_after  dft_after.d2h

export MOL=co2_d2h.mol
export MENU=lda_levyle_nosdft_shielding.inp
#../run.sh
#export MENU=blyp_levyle_nosdft_shielding.inp
#../run.sh
#export MENU=b3lyp_levyle_nosdft_shielding.inp
#../run.sh

export MOL=co2_c1.mol
export MENU=lda_levyle_nosdft_shielding.inp
#../run.sh
export MENU=blyp_levyle_nosdft_shielding.inp
#../run.sh
export MENU=b3lyp_levyle_nosdft_shielding.inp
#../run.sh

export MOL=co2_d2h.mol
export MENU=lda_lvcorr_sdft_shielding.inp
../run.sh
export MENU=lda_lvcorr_nosdft_shielding.inp
../run.sh

#
#for molecule in co2*.mol; do
#    export MOL=${molecule%\.mol}
#
#    export MENU=dhf_lvcorr_shielding.inp
#    #../run.sh
#
#    export MENU=dhf_levyle_shielding.inp
#    #../run.sh
#
#    export MENU=dhf_lvcorr_bzlao.inp
##    ../run.sh
#    export MENU=dhf_levyle_bzlao.inp
##    ../run.sh
#    export MENU=lda_lvcorr_nosdft_shielding.inp
#    ../run.sh
#
#    export MENU=lda_levyle_nosdft_shielding.inp
#    ../run.sh
#
#    export MENU=blyp_levyle_nosdft_shielding.inp
#    ../run.sh
#    export MENU=b3lyp_levyle_nosdft_shielding.inp
#    ../run.sh
#
#    export MENU=lda_lvcorr_sdft_shielding.inp
#    ../run.sh
#
#    export MENU=dhf_lvcorr_magnetizability.inp
##    ../run.sh
#    export MENU=dhf_levyle_magnetizability.inp
##    ../run.sh
#done
#

#!!!!
    export MOL=co2_c1.mol

    export MENU=dhf_lvcorr_shielding.inp
#    ../run.sh
    export MENU=dhf_levyle_shielding.inp
#    ../run.sh
    export MENU=dhf_lvcorr_bzlao.inp
#    ../run.sh
    export MENU=dhf_levyle_bzlao.inp
#    ../run.sh
    export MENU=lda_lvcorr_nosdft_shielding.inp
#    ../run.sh
    export MENU=lda_levyle_nosdft_shielding.inp
#    ../run.sh
    export MENU=blyp_levyle_nosdft_shielding.inp
#    ../run.sh
    export MENU=b3lyp_levyle_nosdft_shielding.inp
#    ../run.sh

    export MENU=dhf_lvcorr_magnetizability.inp
#    ../run.sh
    export MENU=dhf_levyle_magnetizability.inp
#    ../run.sh

exit 0
