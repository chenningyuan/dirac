

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
* improved X2C (reduced memory cost) + local spin-free and spin-orbit X2C [DLU scheme of Peng and Reiher] (C1 symmetry; other double groups perhaps)
* atomic spin-orbit code (spin-same orbit) replacing AMFI
* commutator-based CC (to be discussed with Timo and Mickael)
* PCM (check with Roberto)
* Pol embedding (to see with Erik)


New defaults
------------

* if linked to MKL and nr of MKL threads not set, we quit with clear message what to do
* OMP threads should be set to 1
