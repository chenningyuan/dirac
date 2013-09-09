#!/bin/sh

# set up proper position of your pam launcher
PAM=../../../build_intel_mkl_int64_debug/pam

# F fragment:
$PAM --noarch --inp=F.x2c.scf.inp --mol=F.sto2g.lsym.mol --get "DFACMO=DFFAT"
# Si fragment:
$PAM --noarch --inp=Si.x2c.scf.inp --mol=Si.sto2g.lsym.mol --get "DFACMO=DFSiAT"
# full SiF4  molecule
$PAM --noarch --inp=SiF4.x2c.scf_projanal.inp --mol=SiF4.sto2g.lsym.mol --copy "DFSiAT DFFAT=DFFAT1 DFFAT=DFFAT2 DFFAT=DFFAT3 DFFAT=DFFAT4"

# remove files with MO-coeff
/bin/rm DFFAT DFSiAT

exit
