#!/bin/sh


for molecule in co2*.mol; do
    export MOL=${molecule%\.mol}

#  test dft:
#  ---------

    export MENU=lda_levyle_nosdft_shielding.inp
    ../SCRIPTS/run.sh
    export MENU=blyp_levyle_nosdft_shielding.inp
    ../SCRIPTS/run.sh
    export MENU=b3lyp_levyle_nosdft_shielding.inp
    ../SCRIPTS/run.sh

    export MENU=lda_lvcorr_nosdft_shielding.inp
    #../SCRIPTS/run.sh


#  test sdft:
#  ----------

    export MENU=lda_lvcorr_sdft_shielding.inp
    #../SCRIPTS/run.sh

#   activate to test also dhf:
#   --------------------------
    export MENU=dhf_lvcorr_shielding.inp
    #../SCRIPTS/run.sh

    export MENU=dhf_levyle_shielding.inp
    #../SCRIPTS/run.sh

done

exit 0
