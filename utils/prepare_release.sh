#!/bin/bash

echo "** Dirac08 release preparation script **"
echo

# First make sure that we run from the Dirac root, and that
# the user understands that files will be modified _in place_

if [[ ! -d dirac || ! -d abacus ]]; then
    echo "This script needs to be run from the Dirac source code root."
    exit -1
fi

if [[ -r ../Dirac08 ]]; then
    echo "Directory ../Dirac08 already exist, I cannot proceed."
    exit 1
fi

# Make a copy of the current source tree with tar (to preserve
# symlinks and other properties)
echo "Copying source tree to ../Dirac08"

mkdir ../Dirac08
tar c . | (cd ../Dirac08; tar x)

cd ../Dirac08

#
# First delete code that should not be released.
#
delfunc="-UMOD_ECP -UMOD_KRMC -UMOD_DNF -UMOD_VISUAL -UMOD_CDFT -UMOD_OOF -UMOD_CAP -UMOD_HSCC -UMOD_IHFSCC -UMOD_AMFI -UMOD_LONDON -UMOD_FMM -UMOD_EMB -UMOD_XML"

echo "Removing functionalities from source files"
utils/cppp.py $delfunc `find . -name '*.F' -o -name '*.f90' -o -name '*.F90' -o -name '*.c'`

# Then delete files that should not be shipped
echo "Deleting files"
rm utils/cppp.py utils/prepare_release.sh
rm pamforay Programmers_manual.txt dumbconfigure TODO
rm utils/numerov.py
rm utils/numerov_parameters.py
rm utils/add_copyright.py
rm utils/basis_set.py
rm utils/codegen_maxima.py
rm utils/codegen_mupac.py
rm utils/ihfscc_input_helper.py
rm utils/make_targetbuilds.py
rm utils/periodic_table.py
rm utils/prepare_dyall_basis.sh
rm utils/prepare_dyall_dblock_basis.sh
rm utils/prepare_dyall_fblock_basis.sh
rm utils/prepare_dyall_pblock_basis.sh
rm utils/prepare_dyall_sblock_basis.sh
rm utils/reformat_includes.py
rm utils/remove_copyright.py
rm utils/qjaco.F
rm -rf utils/foray_0.1
rm -rf utils/visualization_scripts
rm -rf notes

# delete tests that are not shipped
# and adapt runtest.in
echo "Deleting test directories that are not shipped"
for test in 43.plots \
            56.aoresponse \
            45.densfit \
            26.krmc_gosci \
            27.krmc_luciarel \
            28.luciarel_medium \
            47.fde.prop \
            38.london.magn_prop; do

    rm -rf test/$test
    sed "s/$test//g" test/runtest.in > delete_me; mv delete_me test/runtest.in
done

for directory in hsmrcc \
                 cap \
                 oof \
                 amfi \
                 argosecp \
                 arducca \
                 visual \
                 cdft \
                 krmc \
                 luciarel \
                 aoresponse \
                 amfi; do

    echo "Deleting directory $directory"
    rm -rf $directory
done

# Delete .svn files if there are any (should not be with svn export)
echo "Deleting old .svn files"
rm -rf `find . -name '.svn'`

#
# Figure out the dependencies that will be shipped.
# This requires that you have configured the local Dirac
# correctly.
#
if [ -r Makefile.config ]; then
    echo "Creating dependencies for release version"
    echo "I will assume you have configured Dirac ok"
    make depend > /dev/null
else
    echo "Cannot find Makefile.config"
    echo "  - Dirac08.tar.gz will NOT include dependencies"
fi


# run pamadm with the correct flags
#
echo "Running pamadm -dft -xpp -eri -adc export"
./pamadm -dft -xpp -eri -adc export > /dev/null
if [[ -r Dirac08.tar.gz ]]; then
    echo "Created ../Dirac08/Dirac08.tar.gz"
else
    echo "Failed to create ../Dirac08/Dirac08.tar.gz"
    echo "Run pamadm in ../Dirac08 to debug"
fi

