#!/bin/sh

#
# $Id$
#
# shell script to prepare intermediate files for assembling dyall basis dirac/dalton basis set files
#

ARCHIVE_DIR=../basis/reference_unedited


zeta=dz

block=4s
val="correlating,4s4p polarizing"
core="correlating,4s4p polarizing"
python2.4 basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                       --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                       --polarization_type="$val" \
                       --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.val.bas  
python2.4 basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                       --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                       --polarization_type="$core" \
                       --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.core.bas  


block=5s
val="correlating,4s4p polarizing"
core="correlating,4s4p polarizing"
python2.4 basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                       --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                       --polarization_type="$val" \
                       --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.val.bas  
python2.4 basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                       --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                       --polarization_type="$core" \
                       --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.core.bas  


block=6s
val="correlating,5s5p polarizing"
core="correlating,5s5p polarizing"
python2.4 basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                       --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                       --polarization_type="$val" \
                       --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.val.bas  
python2.4 basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                       --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                       --polarization_type="$core" \
                       --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.core.bas  


block=7s
val="correlating,4s4p polarizing"
core="correlating,4s4p polarizing"
python2.4 basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                       --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                       --polarization_type="$val" \
                       --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.val.bas  
python2.4 basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                       --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                       --polarization_type="$core" \
                       --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.core.bas  

#
block=4s
val="4s4p correlating,3s3p correlating,3s3p polarizing"
core="4s4p correlating,3s3p correlating,3s3p polarizing,2s2p correlating"
for zeta in tz qz; do
python2.4 basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                       --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                       --polarization_type="$val" \
                       --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.val.bas
python2.4 basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                       --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                       --polarization_type="$core" \
                       --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.core.bas
done

block=5s
val="5s5p correlating,4s4p correlating,4s4p polarizing"
core="5s5p correlating,4s4p correlating,4s4p polarizing,3d correlating"
for zeta in tz qz; do
python2.4 basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                       --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                       --polarization_type="$val" \
                       --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.val.bas  
python2.4 basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                       --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                       --polarization_type="$core" \
                       --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.core.bas  
done

block=6s
val="6s6p correlating,5s5p correlating,5s5p polarizing"
core="6s6p correlating,5s5p correlating,5s5p polarizing,4d correlating"
for zeta in tz qz; do
python2.4 basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                       --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                       --polarization_type="$val" \
                       --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.val.bas  
python2.4 basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                       --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                       --polarization_type="$core" \
                       --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.core.bas  
done

block=7s
val="7s7p correlating,6s6p correlating,6s6p polarizing"
core="7s7p correlating,6s6p correlating,6s6p polarizing,5d correlating"
zeta=tz
python2.4 basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                       --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                       --polarization_type="$val" \
                       --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.val.bas  
python2.4 basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                       --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                       --polarization_type="$core" \
                       --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.core.bas  

block=7s
val="7s7p correlating,6s6p correlating,6s6p dipole polarization"
core="7s7p correlating,6s6p correlating,6s6p dipole polarization,5d correlating"
zeta=qz
python2.4 basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                       --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                       --polarization_type="$val" \
                       --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.val.bas  
python2.4 basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                       --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                       --polarization_type="$core" \
                       --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.core.bas  
