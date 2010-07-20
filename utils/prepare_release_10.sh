#!/bin/bash

source_directory=$PWD

release='DIRAC10'
echo "** $release release preparation script **"
echo


# first make sure that we run from the dirac root, and that
# the user understands that files will be modified _in place_
# ===========================================================

if [[ ! -d dirac || ! -d abacus ]]; then
    echo "This script needs to be run from the Dirac source code root."
    exit -1
fi
if [[ -r ../$release ]]; then
    echo "Directory ../$release already exist, I cannot proceed."
    exit 1
fi


# make a copy of the current source tree with tar (to preserve symlinks and other properties)
# ===========================================================================================

echo "Copying source tree to ../$release"
mkdir ../$release
tar c . | (cd ../$release; tar x)
cd ../$release


# remove files and directories
# ============================

# do not remove cppp.py here - it is needed and removed later
echo "Deleting files that should not be shipped"
for file in \
    utils/prepare_release_08.sh         \
    utils/prepare_release_10.sh         \
    utils/update_copyright.py           \
    doc/release_statement               \
    utils/add_copyright.py              \
    utils/basis_set.py                  \
    utils/codegen_maxima.py             \
    utils/codegen_mupac.py              \
    utils/ihfscc_input_helper.py        \
    utils/make_targetbuilds.py          \
    utils/periodic_table.py             \
    relccsd/ccints_nomemget.F           \
    relccsd/ccdriv_nomemget.F           \
    utils/prepare_dyall_basis.sh        \
    utils/prepare_dyall_dblock_basis.sh \
    utils/prepare_dyall_fblock_basis.sh \
    utils/prepare_dyall_pblock_basis.sh \
    utils/prepare_dyall_sblock_basis.sh \
    utils/reformat_includes.py          \
    utils/qjaco.F                       \
  ; do
    if [ -e "$file" ]; then
       rm -f $file
    else
       echo "WARNING: file $file does not exist"
    fi
done

echo "Deleting directories that should not be shipped"
for directory in \
    utils/foray_0.1             \
    utils/visualization_scripts \
    utils/diracinput            \
    utils/fortran               \
    utils/radovan               \
    notes                       \
    hsmrcc                      \
    cap                         \
    oof                         \
    xcfun                       \
    argosecp                    \
    arducca                     \
    krmc                        \
    luciarel                    \
    openrsp                     \
    prp/esr                     \
    amfi                        \
  ; do
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


# we need some files that have been removed further above
# =======================================================

echo "Copying some necessary files that have been removed further above"
mkdir openrsp
cp $source_directory/openrsp/dirac_interface.F90 openrsp
cp $source_directory/openrsp/precision.F90       openrsp


# remove unreleased functionality (do this after files have been removed or restored)
# ===================================================================================

echo "Removing functionalities from source files that should not be released"
delfunc="-UMOD_UNRELEASED -UMOD_OPENRSP -UMOD_ESR -UMOD_ECP -UMOD_KRMC -UMOD_DNF -UMOD_XCFUN -UMOD_OOF -UMOD_CAP -UMOD_HSCC -UMOD_IHFSCC -UMOD_AMFI -UMOD_FMM -UMOD_EMB -UMOD_XML"
utils/cppp.py $delfunc `find . -name '*.F' -o -name '*.f90' -o -name '*.F90' -o -name '*.c'`
# now remove utils/cppp.py
rm -f utils/cppp.py


# remove and modify some tests
# ============================

echo "Deleting tests that should not be shipped"
for test in \
    openrsp_cars               \
    openrsp_cme                \
    openrsp_delta              \
    openrsp_electric           \
    openrsp_molecular_gradient \
    openrsp_roa                \
  ; do
    if [ -d "test/$test" ]; then
       rm -rf test/$test
       sed "s/$test//g" test/runtest > delete_me; mv delete_me test/runtest
    else
       echo "WARNING: test $test does not exist"
    fi
done

echo "Removing certain functionalities from tests that are shipped"
for test in \
    19.fockspace      \
    6.cosci.energy    \
    10.efg            \
    11.methods.cosci  \
    36.mp2_properties \
    2.nmrspinspin     \
    5.cisd.energy     \
    33.frozen         \
  ; do
    sed '/REMOVE_BEGIN/,/REMOVE_END/'d test/$test/run.sh > test/$test/run.sh_red
    mv test/$test/run.sh_red  test/$test/run.sh
    chmod u+x test/$test/run.sh
done


# figure out the dependencies that will be shipped
# this requires that you have configured the local dirac correctly
# ================================================================

if [ -r Makefile.config ]; then
    echo "Creating dependencies for release version"
    echo "I will assume you have configured Dirac (Makefile.config)"
    make depend > /dev/null
else
    echo "Cannot find Makefile.config"
    echo "  - $release.tar.gz will NOT include dependencies"
fi


# run pamadm with the correct flags
# =================================

echo "Running ./pamadm -xpp -eri -adc -visual -london export"
./pamadm -xpp -eri -adc -visual -london export > /dev/null
if [[ -r $release.tar.gz ]]; then
    echo "Created ../$release/$release.tar.gz"
else
    echo "Failed to create ../$release/$release.tar.gz"
    echo "Run pamadm in ../$release to debug"
fi

