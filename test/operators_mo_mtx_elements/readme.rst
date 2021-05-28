===========================
Operator MO-matrix elements
===========================

Testing few property operators available in DIRAC.

Property MO-matrix elements over occupied orbitals of small relativistic atom are printed out at the post-SCF level 
thanks to the PRPTRA module.

At the two-component (X2c) level, picture change transformed MO elements printing is possible.
Another possible aspect is to test MO-transformed operators at the spin-free level, both 4-component and 2-component.

You launch your operator test as: ::

 ./pam --noarch --inp=Ne.dc_rkb.2fs.scf_prptra_YOUROPERATOR.inp  --mol=Rn_Ne-like.mol

Input files of this test are supposed to be part of the web-tutorial. However, only few selected cases are run here.

Note that MO-matrix elements are checked  only over s- and p- atomic shells here.

One can extend this test with more property operators, which can be contructed in DIRAC.
