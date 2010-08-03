#!/bin/sh

# gosia: this test is the supplement to 38.london.magn_prop
# and to 61.london_dft_magn_prop
# and may be copied there (or removed/renamed...)
# the purpose is to test symmetry in prp/*london* magnetic properties
# with:
#  dhf: lvcorr and levy-leblond 
#  dft:          - || -
#  sdft:         - || -


for molecule in co2*.mol; do
    export MOL=${molecule%\.mol}

    export MENU=dhf_lvcorr_shielding.inp
    ../run.sh
    export MENU=dhf_levyle_shielding.inp
    ../run.sh
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
done


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
