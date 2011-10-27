

New features in DIRAC11
-----------------------

* Analytic molecular gradient at the DFT level
* New and fast XC integration
* Functional derivatives using automatic differentiation (XCFun)
* New visualization options


New build system and infrastructure
-----------------------------------

* New compilation scheme: configure replaced by CMake mechanism
* New pam script (python)
* Alternative launcher: wrapper.py (python)
* New testing framework based on python (runscript)


Important input changes
-----------------------

* XC GRID has own input section
* .DHF is now .SCF


Changed defaults
----------------


Methods
-------

* Hartree-Fock
* Density Functional Theory
* Coupled Cluster
* Configuration Interaction
* Moeller-Plesset Perturbation Theory


Hamiltonians
------------

* 4c Dirac-Coulomb (includes scalar relativistic effects and spin-own-orbit coupling)
* 4c Dirac-Coulomb-Gaunt (includes also spin-other-orbit coupling) (only HF)
* 4c spin-free Dirac-Coulomb (scalar relativistic effects only)
* 4c Levy-Leblond (nonrelativistic)
* 2c X2C, the one-step exact two-component Hamiltonian
* 2c BSS, the two-step exact two-component Hamiltonian (= DKH(infinity,0))


Molecular properties
--------------------

* Up to quadratic response properties at the Hartree-Fock and DFT level
* First order properties with MP2
* Core excitation energies in the static exchange (STEX) approximation.
* Ionization energies at the ADC(3) level of theory


Efficiency
----------

* Full symmetry handling for linear molecules (otherwise up to D2h)
* Parallelization using MPI library calls (MPI should be pre-installed)
