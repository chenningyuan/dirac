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

# First delete code that should not be released.
echo "Removing functionality from source files"
utils/cppp.py -UMOD_KRMC -UMOD_DNF `find . -name '*.F' -o -name '*.f90' -o -name '*.F90' -o -name '*.c'`

# Then delete files that should not be shipped
echo "Deleting files"
rm utils/cppp.py utils/prepare_release.sh
rm pamforay Programmers_manual.txt alternative_pamadm.py dumbconfigure TODO
rm -rf utils/foray_0.1


# Delete .svn files if there are any (should not be with svn export)
echo "Deleting old .svn files"
rm -rf `find . -name '.svn'`

# run pamadm with the correct flags
# 
echo "Running pamadm -dft -xpp export"
./pamadm -dft -xpp export > /dev/null
if [[ -r Dirac08.tar.gz ]]; then
    echo "Created ../Dirac08/Dirac08.tar.gz"
else
    echo "Failed to create ../Dirac08/Dirac08.tar.gz"
    echo "Run pamadm in ../Dirac08 to debug"
fi

