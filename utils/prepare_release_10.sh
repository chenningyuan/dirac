#!/bin/bash

release='DIRAC10'

echo "** $release release preparation script **"
echo

# First make sure that we run from the Dirac root, and that
# the user understands that files will be modified _in place_

if [[ ! -d dirac || ! -d abacus ]]; then
    echo "This script needs to be run from the Dirac source code root."
    exit -1
fi

if [[ -r ../$release ]]; then
    echo "Directory ../$release already exist, I cannot proceed."
    exit 1
fi

# Make a copy of the current source tree with tar (to preserve
# symlinks and other properties)
echo "Copying source tree to ../$release"

mkdir ../$release
tar c . | (cd ../$release; tar x)

cd ../$release

#
# First delete code that should not be released.
#
delfunc="-UMOD_ECP -UMOD_KRMC -UMOD_DNF -UMOD_XCFUN -UMOD_OOF -UMOD_CAP -UMOD_HSCC -UMOD_IHFSCC -UMOD_AMFI -UMOD_LONDON -UMOD_FMM -UMOD_EMB -UMOD_XML"

echo "Removing functionalities from source files"
utils/cppp.py $delfunc `find . -name '*.F' -o -name '*.f90' -o -name '*.F90' -o -name '*.c'`

# Then delete files that should not be shipped
echo "Deleting files"
rm utils/cppp.py utils/prepare_release_10.sh
rm doc/release_statement
rm utils/add_copyright.py
rm utils/basis_set.py
rm utils/codegen_maxima.py
rm utils/codegen_mupac.py
rm utils/ihfscc_input_helper.py
rm utils/make_targetbuilds.py
rm utils/periodic_table.py
rm relccsd/ccints_nomemget.F
rm relccsd/ccdriv_nomemget.F
rm utils/prepare_dyall_basis.sh
rm utils/prepare_dyall_dblock_basis.sh
rm utils/prepare_dyall_fblock_basis.sh
rm utils/prepare_dyall_pblock_basis.sh
rm utils/prepare_dyall_sblock_basis.sh
rm utils/reformat_includes.py
rm utils/remove_copyright.py
rm utils/copyright_fortran
rm utils/copyright_c
rm utils/qjaco.F
rm -rf utils/foray_0.1
rm -rf utils/visualization_scripts
rm -rf utils/diracinput
rm -rf utils/fortran
rm -rf notes
# Delete reference files that can be fetched from the Dirac website
rm basis/reference_unedited/*[^sh]



# delete tests that are not shipped
# and adapt runtest
echo "Deleting test directories that are not shipped"
for test in openrsp_cars \
            openrsp_cme \
            openrsp_delta \
            openrsp_electric \
            openrsp_molecular_gradient \
            openrsp_roa; do

    rm -rf test/$test
    sed "s/$test//g" test/runtest > delete_me; mv delete_me test/runtest
done

echo "Removing certain functionalities from tests"
for this_test in  19.fockspace \
              6.cosci.energy  \
             10.efg  \
             11.methods.cosci  \
             36.mp2_properties  \
              2.nmrspinspin  \
              5.cisd.energy  \
              33.frozen  \
	       ; do

    sed '/REMOVE_BEGIN/,/REMOVE_END/'d test/$this_test/run.sh > test/$this_test/run.sh_red 
    mv test/$this_test/run.sh_red  test/$this_test/run.sh
    chmod u+x test/$this_test/run.sh

done

for directory in hsmrcc \
                 cap \
                 oof \
                 xcfun \
                 argosecp \
                 arducca \
                 krmc \
                 luciarel \
                 openrsp \
                 amfi; do
    echo "Deleting directory $directory"
    rm -rf $directory
done

# Delete junk
echo "Deleting *~ *.o"
rm -f `find . -name '*~' -o -name '*.o'`

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
    echo "I will assume you have configured Dirac (Makefile.config)"
    make depend > /dev/null
else
    echo "Cannot find Makefile.config"
    echo "  - $release.tar.gz will NOT include dependencies"
fi


# run pamadm with the correct flags
#
echo "Running pamadm -xpp -eri -adc -visual export"
./pamadm -xpp -eri -adc -visual export > /dev/null
if [[ -r $release.tar.gz ]]; then
    echo "Created ../$release/$release.tar.gz"
else
    echo "Failed to create ../$release/$release.tar.gz"
    echo "Run pamadm in ../$release to debug"
fi

