#!/bin/sh

#
# template for how to prepare basis set library files for the dyall basis
#
# andre gomes, 2006-10-06
#
# 1. first we check what kind of polarization functions there are in the basis set files

ARCHIVE_DIR=../basis/reference_unedited

#
# 2. taking as an example the 'valence correlating' type of functions, we create the files
#

for block in  4p 5p 6p; do

   zeta=dz

   python2.4 basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive.txt" \
                          --print_polarization --print_scfbas --uncontracted \
                          --polarization_type='correlating' \
                          --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.val.bas  

   python2.4 basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive.txt" \
                          --print_polarization --print_scfbas --uncontracted \
                          --polarization_type='correlating' \
                          --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.core.bas  


   python2.4 basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive.txt" \
                          --print_polarization --print_scfbas --uncontracted \
                          --polarization_type='correlating','diffuse' \
                          --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.aug-val.bas  

   python2.4 basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive.txt" \
                          --print_polarization --print_scfbas --uncontracted \
                          --polarization_type='correlating','diffuse' \
                          --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.aug-core.bas  

   for zeta in tz qz; do

      python2.4 basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive.txt" \
                             --print_polarization --print_scfbas --uncontracted \
                             --polarization_type='valence correlating' \
                             --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.val.bas  


      python2.4 basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive.txt" \
                             --print_polarization --print_scfbas --uncontracted \
                             --polarization_type='valence correlating','core correlating' \
                             --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.core.bas  

      python2.4 basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive.txt" \
                             --print_polarization --print_scfbas --uncontracted \
                             --polarization_type='valence correlating','diffuse' \
                             --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.aug-val.bas  


      python2.4 basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive.txt" \
                             --print_polarization --print_scfbas --uncontracted \
                             --polarization_type='valence correlating','core correlating','diffuse' \
                             --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.aug-corel.bas  
   done
done

