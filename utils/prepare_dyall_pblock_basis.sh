#!/bin/sh

#
# $Id$
#
# shell script to prepare intermediate files for assembling dyall basis dirac/dalton basis set files
#

ARCHIVE_DIR=../basis/reference_unedited

for block in  4p 5p 6p; do

#note, for the double-zeta, add the polarizing f function by hand on the core.bas files?
   zeta=dz

   python2.4 basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                          --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                          --polarization_type='correlating' \
                          --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.val.bas  

   python2.4 basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                          --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                          --polarization_type='correlating' \
                          --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.core.bas  


   python2.4 basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                          --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                          --polarization_type='correlating','diffuse' \
                          --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.aug-val.bas  

   python2.4 basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                          --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                          --polarization_type='correlating','diffuse' \
                          --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.aug-core.bas  

   for zeta in tz qz; do

      python2.4 basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                             --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                             --polarization_type='valence correlating' \
                             --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.val.bas  


      python2.4 basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                             --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                             --polarization_type='valence correlating','core correlating' \
                             --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.core.bas  

      python2.4 basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                             --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                             --polarization_type='valence correlating','diffuse' \
                             --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.aug-val.bas  


      python2.4 basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                             --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                             --polarization_type='valence correlating','core correlating','diffuse' \
                             --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.aug-core.bas  
   done
done

