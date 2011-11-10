

New features in DIRAC11
-----------------------

* Analytic molecular gradient at the DFT level
* New and fast XC integration
* Functional derivatives using automatic differentiation (XCFun)
* New visualization options
* RKBIMP: MO-coefficients generated using restricted kinetic balance (RKB) can
  be extended by their unrestricted kinetic balance (UKB) complement, thus
  providing magnetic balance for response calculations involving external
  magnetic fields
* New and improved 2c Hamiltonian schemes


New build system and infrastructure
-----------------------------------

* New compilation scheme: configure replaced by CMake mechanism
* New pam script (python)
* Alternative launcher: wrapper.py (python)
* New testing framework based on python (runscript)
* Many static allocation calls replaced by dynamic allocation; in practice
  this means that you may need less WORK array memory and/or more space for
  dynamic allocation compared to DIRAC10.


Important input changes
-----------------------

* XC GRID has own input section
* .DHF is now .SCF


Changed defaults
----------------

* .LVCORR is now default
  you can force explicit evaluation of (SS|SS) integrals with .DOSSSS


Methods
-------

* Hartree-Fock
* Density Functional Theory
* Kramers-restricted Multi-Configuration Self-Consistent-Field
* Coupled Cluster
* Configuration Interaction
* Moeller-Plesset Perturbation Theory


Hamiltonians
------------

* 4c Dirac-Coulomb (includes scalar relativistic effects and spin-own-orbit coupling)
* 4c Dirac-Coulomb-Gaunt (only HF; includes also spin-other-orbit coupling)
* 4c spin-free Dirac-Coulomb (scalar relativistic effects only)
* 4c Levy-Leblond (nonrelativistic)
* 2c X2C, the one-step exact two-component Hamiltonian
* 2c BSS, the two-step exact two-component Hamiltonian (= DKH(infinity,0))
* 2c molecular-mean-field (= X2Cmmf), X2C transformation with the converged 4c-Fock operator as defining Hamiltonian


Molecular properties
--------------------

* Up to quadratic response properties at the HF and DFT level
* First-order properties with MP2
* Core excitation energies in the static exchange (STEX) approximation
* Ionization energies at the ADC(3) level of theory
* Selected first-order properties with CI


Efficiency
----------

* Full symmetry handling for linear molecules (otherwise up to D2h)
* Parallelization using MPI library calls (MPI should be pre-installed)
