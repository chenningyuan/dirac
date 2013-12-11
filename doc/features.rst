

New features in DIRAC13
-----------------------

* magnetizabilities at the Hartree-Fock level with London atomic orbitals (LAOs)
* `improved output for TDDFT excitation energies <tutorials/dft/tddft.html>`_
* Hueckel start with atomic fragments
* XML output functionality
* enhancements to frozen density embedding (`FDE <manual/fde.html>`_)
* polarization propagator for 4-C excitations (ADC2 extended)
* enhancements to X2C: local spin-free and spin-orbit X2C (`X2C and local X2C <tutorials/two_component_hamiltonians/x2c_mol_loc.html>`_)
* Dyall basis sets redefined to reduce linear dependence and conform to basis archive files, including fixes
* Basis sets for 2s, 2p, 3s, 3p, 3d added to Dyall 2z, Dyall 3z sets and for 2p, 3s, 3p for Dyall 4z basis sets.
* Dyall aenz (all-electron) basis sets added, with correlating functions for all shells.
* works on Windows 7/8 with GNU MinGW32/64 suite and with native math libraries (see `installation <installation/windows.html>`_)


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
