==========
DIRAC-MRCC
==========

This test is meant to test the interface to the MRCC code of Kallay and coworkers.

It does a standard CCSD calculation with DIRAC - RELCCSD and then produces the interface file fort.55
and an input file fort.56 for the MRCC suite of codes. The two codes should obviously produce the
same CCSD energy, which is the purpose of this test.

NB: Since MRCC is not distributed with DIRAC, you should get your own copy of the code from
www.mrcc.hu and set the environment variable $MRCCDIR to point to the "goldstone" and
"mrcc" executables. With this set, the test script will execute correctly, but it can for 
technical reasons (setup of the test scripts) not yet automatically
test the outcome of the MRCC calculation. Please make sure to check the output yourself.

Example of the run
~~~~~~~~~~~~~~~~~~
export MRCCDIR=/home/milias/work/software/mrcc  # position of MRCC code from mrcc.hu
export DIRAC_UTILS_PATH=/home/milias/work/software/dirac/trunk_cloned/build_gnu_mkl_ilp64 #  dirac.x, dirac_mo_integrls_export.x

./test -b /home/milias/work/software/dirac/trunk_cloned/build_gnu_mkl_ilp64 > test.log
