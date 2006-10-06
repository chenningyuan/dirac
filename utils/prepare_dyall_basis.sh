#!/bin/sh

#
# template for how to prepare basis set library files for the dyall basis
#
# andre gomes, 2006-10-06
#
# 1. first we check what kind of polarization functions there are in the basis set files

ARCHIVE_DIR=../basis/reference_unedited

for block in  4p 4d 5s 5p 5d 5f 6p; do
   for zeta in dz tz qz; do

      python2.4 basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive.txt" \
                             --list_polarization > $block"_"$zeta"_polarization-types.txt"

   done
done

#
# 2. taking as an example the 'valence correlating' type of functions, we create the files
#

for block in  4p 4d 5s 5p 5d 5f 6p; do
   for zeta in tz qz; do

      python2.4 basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive.txt" \
                             --print_polarization --print_scfbas --uncontracted \
                             --polarization_type='valence correlating' \
                             --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.bas  

   done
done

#
# 3. concatenate the resulting intermediate .bas files by hand, and edit the "elements supported field"
#

for zeta in tz qz; do
   cat dyall*.$zeta.bas > dyall.valence_correlating.$zeta.bas
done 

#
# 4. repeat 2 and 3 above for other types of correlating functions
#
