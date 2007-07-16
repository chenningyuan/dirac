#!/bin/sh

#
# $Id$
#
# shell script to prepare intermediate files for assembling dyall basis dirac/dalton basis set files
#

ARCHIVE_DIR=../basis/reference_unedited

block=5f
zeta=dz

val="5f correlating,6d/7s correlating,5f polarizing"
core1="5f correlating,6d/7s correlating,5d correlating,5f polarizing"
core2="5f correlating,6d/7s correlating,5d correlating,5f polarizing,6d polarizing"

python2.4 basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                       --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                       --polarization_type="$val" \
                       --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.val.bas  


python2.4 basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                       --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                       --polarization_type="$core1" \
                       --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.core1.bas  


python2.4 basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                       --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                       --polarization_type="$core2" \
                       --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.core2.bas  

block=5f
zeta=tz

val='5f correlating,6d/7s correlating g,5f polarizing'
core1="5f correlating,6d/7s correlating g,5d correlating,5f polarizing"  
core2="5f correlating,6d/7s correlating g,5d correlating,5f polarizing,6d polarizing"

echo $val
python2.4 basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                       --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                       --polarization_type="$val" \
                       --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.val.bas  


python2.4 basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                       --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                       --polarization_type="$core1" \
                       --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.core1.bas  


python2.4 basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                       --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                       --polarization_type="$core2" \
                       --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.core2.bas  

block=5f
zeta=qz


val="5f correlating,6d/7s correlating,6s6p correlating,5f polarizing"
core1="5f correlating,6d/7s correlating,6s6p correlating,5d correlating,5f polarizing"
core2="5f correlating,6d/7s correlating,6s6p correlating,5d correlating,5f polarizing,6d/7s polarizing"

python2.4 basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                       --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                       --polarization_type="$val" \
                       --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.val.bas  


python2.4 basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                       --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                       --polarization_type="$core1" \
                       --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.core1.bas  


python2.4 basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                       --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                       --polarization_type="$core2" \
                       --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.core2.bas  

