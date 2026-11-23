:orphan:
 

CFREAD
======

CFREAD is a utility program for reading the legacy binary DFCOEF file that was used
in DIRAC versions prior to DIRAC22. Since DIRAC22, the default checkpoint file is
``CHECKPOINT.h5`` (HDF5 format), which is machine-independent and does not require CFREAD.
CFREAD should therefore be considered a legacy tool for working with old DFCOEF files.

For modern CHECKPOINT.h5 files, use standard HDF5 tools or the appropriate DIRAC analysis
modules instead.

For reference, CFREAD usage on an old DFCOEF file is as follows::

  dirac 26>cfread.x DFCOEF
  
   CFREAD; Information read from file : DFCOEF  (legacy pre-DIRAC22 format)
  
   - Heading :DIRAC: No title specified !!!                     Tue Sep 24 23:05:49 2013
   - Total energy: -75.1621606438559979
   - Fermion ircops          :      1
  -----------------------------------
   - Negative-energy orbitals:      0
   - Positive-energy orbitals:     21
   - AO-basis functions      :     21

Note that the electronic energy is calculated before diagonalization, so that coefficients from iteration *n* 
refer to electronic energy from iteration *n-1*. At convergence this should not matter.

