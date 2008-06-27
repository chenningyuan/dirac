#!/bin/bash

echo "** Dirac08 release preparation script **"
echo

# First make sure that we run from the Dirac root, and that
# the user understands that files will be modified _in place_

if [[ ! -d dirac || ! -d abacus ]]; then
    echo "This script needs to be run from the Dirac source code root."
    exit -1
fi

echo "This will delete code and files in place! Make sure you run this script"
echo "on a copy of your source tree. Please write the word 'continue' to proceed"
read cont
if [[ "$cont" != "continue" ]]; then
    echo "Aborting."
    exit -1
fi

# First delete code that should not be released.
echo "Removing functionality from source files"
utils/cppp.py -UMOD_KRMC -UMOD_DNF `find . -name '*.F' -o -name '*.f90' -o -name '*.F90' -o -name '*.c'`

# Then delete files that should not be shipped
echo "Deleting files"
rm utils/cppp.py utils/prepare_release.sh

# Delete .svn files if there are any (should not be with svn export)
echo "Deleting old .svn files"
rm -rf `find . -name '.svn'`

