

New features in DIRAC13
-----------------------

* magnetizabilities LAO, HF (perhaps DFT)
* improved output for TDDFT excitation energies
* projection analysis in molecular sym (perhaps)
* Hueckel start with atomic fragments
* XML output functionality
* localization of virtuals (Stano, to consider)
* enhancements to FDE
* polarization propagator for 4-C excitations (ADC2 extended)
* improved X2C: local spin-free and spin-orbit X2C [DLU scheme of Peng and Reiher] (C1 symmetry only)
* commutator-based CC (to be discussed with Timo and Mickael)
* Pol embedding (to see with Erik)
* Dyall basis sets redefined to reduce linear dependence and conform to basis archive files, including fixes
* Basis sets for 2s, 2p, 3s, 3p, 3d added to Dyall 2z, Dyall 3z sets and for 2p, 3s, 3p for Dyall 4z basis sets.
* Dyall aenz (all-electron) basis sets added, with correlating functions for all shells.


Improved compilation and testing
--------------------------------

* New test script. Simplified testing using MPI.
* Updated math library detection.
* Better support for MKL libraries.
* Support for Cray.
* Support for MPI runs which do not use mpirun.


New defaults
------------

* The pam script sets MKL_NUM_THREADS=1, MKL_DYNAMIC="FALSE",
  OMP_NUM_THREADS=1, and OMP_DYNAMIC="FALSE" unless these variables are set by
  the user.
