#!/usr/bin/env bash

# motivation for this script is to create a release tarball with all external
# git submodules included so that we can upload a self-contained archive to
# Zenodo and so that users who download the tarball from Zenodo do not have to
# run additional steps

# script does the following steps:
# - create temporary folder
# - git clone --recursive a specific branch
# - strip away git history
# - create tarball
# - clean up temporary folder

set -euf -o pipefail

# adapt the following two lines
VERSION="26.0"
BRANCH="release-26"

THIS_FOLDER=${PWD}
TEMP_FOLDER=$(mktemp -d)

cd ${TEMP_FOLDER}
git clone --recursive --branch ${BRANCH} https://gitlab.com/dirac/dirac.git DIRAC-${VERSION}-Source
mv DIRAC-${VERSION}-Source/.git git-folder-not-shipped
tar -czf DIRAC-${VERSION}-Source.tar.gz DIRAC-${VERSION}-Source
cp DIRAC-${VERSION}-Source.tar.gz ${THIS_FOLDER}

# for easybuild
tar -cJf DIRAC-${VERSION}-Source.tar.xz DIRAC-${VERSION}-Source  
cp DIRAC-${VERSION}-Source.tar.xz ${THIS_FOLDER}

cd ${THIS_FOLDER}
mv ${TEMP_FOLDER} ${THIS_FOLDER}/temp-folder-can-be-removed
rm -rf temp-folder-can-be-removed
