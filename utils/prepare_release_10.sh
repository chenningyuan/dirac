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

echo "Removing functionalities from source files htat should not be released"
delfunc="-UMOD_UNRELEASED -UMOD_ECP -UMOD_KRMC -UMOD_DNF -UMOD_XCFUN -UMOD_OOF -UMOD_CAP -UMOD_HSCC -UMOD_IHFSCC -UMOD_AMFI -UMOD_LONDON -UMOD_FMM -UMOD_EMB -UMOD_XML"
utils/cppp.py $delfunc `find . -name '*.F' -o -name '*.f90' -o -name '*.F90' -o -name '*.c'`

echo "Deleting files that should not be shipped"
for file in utils/prepare_release_08.sh \
            utils/prepare_release_10.sh \
            utils/update_copyright.py \
            utils/cppp.py \
            doc/release_statement \
            utils/add_copyright.py \
            utils/basis_set.py \
            utils/codegen_maxima.py \
            utils/codegen_mupac.py \
            utils/ihfscc_input_helper.py \
            utils/make_targetbuilds.py \
            utils/periodic_table.py \
            relccsd/ccints_nomemget.F \
            relccsd/ccdriv_nomemget.F \
            utils/prepare_dyall_basis.sh \
            utils/prepare_dyall_dblock_basis.sh \
            utils/prepare_dyall_fblock_basis.sh \
            utils/prepare_dyall_pblock_basis.sh \
            utils/prepare_dyall_sblock_basis.sh \
            utils/reformat_includes.py \
            utils/remove_copyright.py \
            utils/copyright_fortran \
            utils/copyright_c \
            utils/qjaco.F; do
    if [ -e "$file" ]; then
       rm -f $file
    else
       echo "WARNING: file $file does not exist"
    fi
done

echo "Deleting directories that should not be shipped"
for directory in utils/foray_0.1 \
                 utils/visualization_scripts \
                 utils/diracinput \
                 utils/fortran \
                 utils/radovan \
                 notes \
                 hsmrcc \
                 cap \
                 oof \
                 xcfun \
                 argosecp \
                 arducca \
                 krmc \
                 luciarel \
                 openrsp \
                 amfi; do
    if [ -d "$directory" ]; then
       rm -rf $directory
    else
       echo "WARNING: directory $directory does not exist"
    fi
done

echo "Deleting reference files that can be fetched from the Dirac website"
for file in `ls basis/reference_unedited/*[^sh]`; do
    rm -f $file
done

echo "Deleting tests that should not be shipped"
for test in openrsp_cars \
            openrsp_cme \
            openrsp_delta \
            openrsp_electric \
            openrsp_molecular_gradient \
            openrsp_roa; do
    if [ -d "test/$test" ]; then
       rm -rf test/$test
       sed "s/$test//g" test/runtest > delete_me; mv delete_me test/runtest
    else
       echo "WARNING: test $test does not exist"
    fi
done

echo "Deleting object files (if any)"
rm -f `find . -name '*.o'`

echo "Deleting junk files (if any)"
for file in `find . -name '*~'`; do
    echo "deleting file $file"
    rm -f $file
done

echo "Deleting .svn directories (if any)"
for directory in `find . -name '.svn'`; do
    echo "deleting directory $directory"
    rm -rf $directory
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

