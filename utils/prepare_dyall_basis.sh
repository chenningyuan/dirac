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
cat *.{4f,4d,4p,5s,5p,5d,5f,6s,6p,6d,7s}.dz.val.bas > dyall.v2z
cat *.{4f,4d,4p,5s,5p,5d,5f,6s,6p,6d,7s}.tz.val.bas > dyall.v3z
cat *.{4f,4d,4p,5s,5p,5d,5f,6s,6p,6d,7s}.qz.val.bas > dyall.v4z

cat *.{4f,4d,4p,5s,5p,5d,5f,6s,6p,6d,7s}.dz.core.bas > dyall.cv2z
cat *.{4f,4d,4p,5s,5p,5d,5f,6s,6p,6d,7s}.tz.core.bas > dyall.cv3z
cat *.{4f,4d,4p,5s,5p,5d,5f,6s,6p,6d,7s}.qz.core.bas > dyall.cv4z

# one must now edit the dyall.{a{,c}}v{2,3,4}z files and
# set the "Elements supported" section correctly by hand.
# in order to do that, one has to just collect all the 
# separate "Elements supported" fields, and leave just one
# with all the atoms
#

#
# as of the time of the latest commit for this file (see Id above)
# the following can be used as template to add the list of supported
# elements to the dyall.{,c}v{2,3,4}z files above:
#
#<---- starts here
#  $
#  $ Elements Supported
#  $ Sr Rb Cs Ba Fr Ra
#  $ Ga As Kr Ge Br Se I Xe Sn In Sb Te Bi Pb Tl At Rn Po
#  $ Ru Ag Nb Tc Pd Y Cd Rh Mo Zr Pt Hf Ir Re Au W Os Hg Ta
#  $ Hs Db Mt Rf Rg Uub Sg Ds Bh
#  $ Md Ac Cm Pu No Am Cf Bk Pa U Th Np Lr Fm Es
#  $ Pr La Yb Nd Ce Gd Ho Sm Dy Eu Tm Er Tb Lu Pm
#  $
#  $  REFERENCE
#  $
#  $ 5s-7s: K.G. Dyall, in preparation.
#  $ 4p-6p: K.G. Dyall, Theor. Chem. Acc. (1998) 99:366;
#  $        addendum    Theor. Chem. Acc. (2002) 108:365;
#  $        revision    Theor. Chem. Acc. (2006) 115:441
#  $    4d: K.G. Dyall, Theor. Chem. Acc. (2007) 117:483.
#  $    5d: K.G. Dyall, Theor. Chem. Acc. 112, 403-409 (2004)
#  $    6d: K.G. Dyall, in preparation.
#  $    5f: K.G. Dyall, Theor. Chem. Acc. (2007) 117:491.
#  $    4f: A.S.P. Gomes and K.G. Dyall, in preparation.
#  $ Available from the Dirac web site, http://dirac.chem.sdu.dk.
#  $
#  $
#<--- ends here
#
# to use this, remove the first 3 characters ("#  ") so that the
# dollar is at the first column (that is, use e.g in vi ":%s/#  //")
#

