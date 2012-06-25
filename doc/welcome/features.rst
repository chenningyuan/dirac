

Features in DIRAC |version|
===========================


New features in DIRAC12 (released 12/12/12)
-------------------------------------------

-  Improved interface with PyADF (distributed separately) - in particular to allow automatized "atomic start guess" SCF runs (DFT and HF)
-  Overlap diagnostic
-  Simple magnetic balance
-  LAOs at the DFT level
-  LAO current densities
-  Fock-space intensities
-  Localization scheme
-  MP2-sr-DFT
-  Atomic start
-  Hückel start
-  Improved projection analysis
-  MP2 natural orbital module
-  CPP
-  New Lanczos
-  XCInt/XCGrid (Tentative)
-  ECP? (we did not discuss this but we should discuss it)

-  New defaults

   -  New input style for RELCC and RELADC
   -  Changed level shift
   -  Changed bare nucleus corrections (new parameters)
   -  New MPI 64/32-interface
   -  Improved start guess and improved SCF convergence
   -  Changed open-shell energy definition


New features in DIRAC11 (released 11/11/11)
-------------------------------------------

-  Analytic molecular gradient at the DFT level
-  New and fast XC integration
-  Functional derivatives using automatic differentiation (XCFun)
-  New visualization options
-  RKBIMP: MO-coefficients generated using restricted kinetic balance
   (RKB) can be extended by their unrestricted kinetic balance (UKB)
   complement, thus providing magnetic balance for response calculations
   involving external magnetic fields
-  New and improved 2c Hamiltonian schemes

-  New build system and infrastructure

   -  New compilation scheme: configure replaced by CMake mechanism
   -  New pam script (python)
   -  Alternative launcher: wrapper.py (python)
   -  New testing framework based on python (runscript)
   -  Many static allocation calls replaced by dynamic allocation; in
      practice this means that you may need less WORK array memory
      and/or more space for dynamic allocation compared to DIRAC10.

-  New mailing lists for help and announcements about the DIRAC program
   package

   -  `http://groups.google.com/group/dirac-users <http://groups.google.com/group/dirac-users>`_

-  Important input changes

   -  XC GRID has own input section
   -  .DHF is now .SCF

-  Changed defaults

   -  .LVCORR is now default; you can force explicit evaluation of
      (SS\|SS) integrals with .DOSSSS

-  Methods

   -  Hartree-Fock
   -  Density Functional Theory
   -  Kramers-restricted Multi-Configuration Self-Consistent-Field
   -  Coupled Cluster
   -  Configuration Interaction
   -  Moeller-Plesset Perturbation Theory

-  Hamiltonians

   -  4c Dirac-Coulomb (includes scalar relativistic effects and
      spin-own-orbit coupling)
   -  4c Dirac-Coulomb-Gaunt (only HF; includes also spin-other-orbit
      coupling)
   -  4c spin-free Dirac-Coulomb (scalar relativistic effects only)
   -  4c Levy-Leblond (nonrelativistic)
   -  2c X2C, the one-step exact two-component Hamiltonian
   -  2c BSS, the two-step exact two-component Hamiltonian (=
      DKH(infinity,0))
   -  2c molecular-mean-field (= X2Cmmf), X2C transformation with the
      converged 4c-Fock operator as defining Hamiltonian

-  Molecular properties

   -  Up to quadratic response properties at the HF and DFT level
   -  First-order properties with MP2
   -  Core excitation energies in the static exchange (STEX)
      approximation
   -  Ionization energies at the ADC(3) level of theory
   -  Selected first-order properties with CI

-  Efficiency

   -  Full symmetry handling for linear molecules (otherwise up to D2h)
   -  Parallelization using MPI library calls (MPI should be
      pre-installed)


New features in DIRAC10
-----------------------

-  Methods:

   -  Kramers-restricted MCSCF
   -  RELADC for correlated calculations of single/double ionization
      spectra
   -  large-scale parallel CI (LUCITA/KRCI)
   -  intermediate Hamiltonian formalism for Fock-space CCSD
   -  interface to MRCC (see
      `http://tc03.fkt.bme.hu/ <http://tc03.fkt.bme.hu/>`_)
   -  frozen density embedding

-  Hamiltonians:

   -  2c X2C+AMFI for 2-electron spin-orbit corrections (spin-same
      orbit[SSO]/spin other-orbit[SOO])

-  Molecular properties:

   -  HF/KS excitation energies
   -  KS response with noncollinear spin polarization and full
      derivative of functionals
   -  linear response functions at imaginary frequencies
   -  more efficient KS DFT code
   -  London orbitals for HF NMR shieldings

-  Analysis:

   -  visualization of unperturbed and perturbed densities
   -  projection analysis of expectation values
   -  expectation values/transition moments KRCI/GOSCI


Features in DIRAC08
-------------------

-  Methods:

   -  Hartree-Fock
   -  Density Functional Theory
   -  Coupled Cluster
   -  Configuration Interaction
   -  Second order Møller-Plesset Perturbation Theory

-  Hamiltonians:

   -  4c Dirac-Coulomb (includes scalar relativistic effects and
      spin-own-orbit coupling)
   -  4c Dirac-Coulomb-Gaunt (includes also spin-other-orbit coupling)
      (only HF)
   -  4c spin-free Dirac-Coulomb (scalar relativistic effects only)
   -  4c Levy-Leblond (nonrelativistic)
   -  2c X2C, the one-step exact two-component Hamiltonian
   -  2c BSS, the two-step exact two-component Hamiltonian (=
      DKH(infinity,0))

-  Molecular properties:

   -  Up to quadratic response properties at the Hartree-Fock and DFT
      level
   -  First order properties with MP2
   -  Core excitation energies in the static exchange (STEX)
      approximation.
   -  Single/Double Ionization energies and spectra at the
      ADC(3)/ADC(2x) level of theory.

-  Efficiency:

   -  Full symmetry handling for linear molecules (otherwise up to D2h)
   -  Parallelization using MPI library calls (MPI should be
      preinstalled)

-  Some of the new features of DIRAC08:

   -  A one-step exact two-component Hamiltonian (X2C)
   -  Relativistic Green's function (propagator) module RELADC for the
      calculation of ionization energies
   -  Possibility to include the Gaunt interaction in HF calculations
   -  Implementation of several new density functionals
   -  Linear and quadratic response DFT
   -  Addition of the latest Dyall basis sets and more non-relativistic
      basis sets to the basis library
   -  Analysis by means of fragment orbitals
   -  New parallelization of the MOLTRA module with reduced I/O
   -  Parallelization of the LUCITA CI module
