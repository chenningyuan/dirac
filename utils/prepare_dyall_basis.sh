#!/bin/sh

#
# $Id$
#
# top-level script to generate atchives for the dyall sets
#

# preparation for the 5s/6s/7s elements
sh prepare_dyall_sblock_basis.sh

# preparation for the 4p/5p/6p elements
sh prepare_dyall_pblock_basis.sh

# preparation for the 4d/5d elements
sh prepare_dyall_dblock_basis.sh

# preparation for the 5f elements
sh prepare_dyall_fblock_basis.sh

# concatenating different blocks, first for sets with diffuses 
cat *.{4p,5p,6p}.dz.aug-val.bas > dyall.av2z
cat *.{4p,5p,6p}.tz.aug-val.bas > dyall.av3z
cat *.{4p,5p,6p}.qz.aug-val.bas > dyall.av4z

cat *.{4p,5p,6p}.dz.aug-core.bas > dyall.acv2z
cat *.{4p,5p,6p}.tz.aug-core.bas > dyall.acv3z
cat *.{4p,5p,6p}.qz.aug-core.bas > dyall.acv4z

# and now withouth diffuses
cat *.{4d,4p,5s,5p,5d,5f,6s,6p,6d,7s}.dz.val.bas > dyall.v2z
cat *.{4d,4p,5s,5p,5d,5f,6s,6p,6d,7s}.tz.val.bas > dyall.v3z
cat *.{4d,4p,5s,5p,5d,5f,6s,6p,6d,7s}.qz.val.bas > dyall.v4z

cat *.{4d,4p,5s,5p,5d,6s,6p,6d,7s}.dz.core.bas > dyall.cv2z
cat *.{4d,4p,5s,5p,5d,6s,6p,6d,7s}.tz.core.bas > dyall.cv3z
cat *.{4d,4p,5s,5p,5d,6s,6p,6d,7s}.qz.core.bas > dyall.cv4z

cat *.5f.dz.core1.bas > dyall.5f.c1v2z
cat *.5f.tz.core1.bas > dyall.5f.c1v3z
cat *.5f.qz.core1.bas > dyall.5f.c1v4z

cat *.5f.dz.core2.bas > dyall.5f.c2v2z
cat *.5f.tz.core2.bas > dyall.5f.c2v3z
cat *.5f.qz.core2.bas > dyall.5f.c2v4z


# one must now edit the dyall.{a{,c}}v{2,3,4}z files and
# set the "Elements supported" section correctly by hand.
# in order to do that, one has to just collect all the 
# separate "Elements supported" fields, and leave just one
# with all the atoms
#
# in particular, one should be careful with any archive that
# contains 5f elements
