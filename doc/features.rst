

New features in DIRAC13
-----------------------

* magnetizabilities LAO, HF (perhaps DFT)
* improved output for TDDFT excitation energies
* projection analysis in molecular sym (perhaps)
* Hueckel start with atomic fragments
* doublet EPR (perhaps)
* XML output functionality
* localization of virtuals (Stano, to consider)
* correction of overlap diagnostic
* enhancements to FDE
* polarization propagator for 4-C excitations (ADC2 extended)
* Fano ADC (perhaps)
* C1 X2C and new AMFI (Stefan please fill in details)
* commutator-based CC (to be discussed with Timo and Mickael)
* PCM (check with Roberto)
* Pol embedding (to see with Erik)


New defaults
------------

* The pam script sets MKL_NUM_THREADS=1, MKL_DYNAMIC="FALSE",
  OMP_NUM_THREADS=1, and OMP_DYNAMIC="FALSE" unless these variables are set by
  the user.
