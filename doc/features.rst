

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
* improved X2C (reduced memory cost) + local spin-free and spin-orbit X2C [DLU scheme of Peng and Reiher] (C1 symmetry; other double groups perhaps)
* atomic spin-orbit code (spin-same orbit) replacing AMFI
* commutator-based CC (to be discussed with Timo and Mickael)
* PCM (check with Roberto)
* Pol embedding (to see with Erik)
* Dyall basis sets redefined to reduce linear dependence and conform to basis archive files, including fixes
* Basis sets for 1s, 2s, 2p, 3s, 3p, 3d added to Dyall *2z, *3z, and *4z basis sets.
* Dyall ae*z (all-electron)  basis sets added, with correlating functions for all shells.

New defaults
------------

* The pam script sets MKL_NUM_THREADS=1, MKL_DYNAMIC="FALSE",
  OMP_NUM_THREADS=1, and OMP_DYNAMIC="FALSE" unless these variables are set by
  the user.
