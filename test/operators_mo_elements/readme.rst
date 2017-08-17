Operator MO-matrix elements
===========================

Test aimed at various property operators available within DIRAC.
MO-matrix elements over occupied orbitals of Xe atoms are printed out at the post-SCF level 
thanks to the PRPTRA module.

At the two-component (X2c) level, picture change transformed MO elements printing is possible.
Another possible aspect is to test MO-transformed operators at the spin-free level, both 4-component and 2-component.

You launch your operatr test as: ::

 ./pam --noarch --put "DFPCMO.Xe.dcrkb.scf-30-24.dyallv2z.lsym=DFPCMO" --inp=Xe.dc_rkb.2fs.scf_prptra_YOUROPERATOR.inp  --mol=Xe.dyall-v2z.lsym.mol

Input files of this test are supposed to be part of the web-tutorial. For this test, however, only few selected cases are run.
