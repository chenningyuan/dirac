

New features in DIRAC11
=======================

* A one-step exact two-component Hamiltonian (X2C) 
* Relativistic Green's function (propagator) module RELADC for the calculation of ionization energies 
* Possibility to include the Gaunt interaction in HF calculations 
* Implementation of several new density functionals 
* Linear and quadratic response DFT 
* Addition of the latest Dyall basis sets and more non-relativistic basis sets to the basis library 
* Analysis by means of fragment orbitals 
* New parallelization of the MOLTRA module with reduced I/O 
* Parallelization of the LUCITA CI module 


Some of the most important regular features
===========================================

Methods
-------

* Hartree-Fock 
* Density Functional Theory 
* Coupled Cluster 
* Configuration Interaction 
* Moeller-Plesset Perturbation Theory 


Hamiltonians
------------

 4c Dirac-Coulomb (includes scalar relativistic effects and spin-own-orbit coupling) 
 4c Dirac-Coulomb-Gaunt (includes also spin-other-orbit coupling) (only HF) 
 4c spin-free Dirac-Coulomb (scalar relativistic effects only) 
 4c Levy-Leblond (nonrelativistic) 
 2c X2C, the one-step exact two-component Hamiltonian 
 2c BSS, the two-step exact two-component Hamiltonian (= DKH(infinity,0)) 


Molecular properties
--------------------

* Up to quadratic response properties at the Hartree-Fock and DFT level 
* First order properties with MP2 
* Core excitation energies in the static exchange (STEX) approximation. 
* Ionization energies at the ADC(3) level of theory. 


Efficiency
----------

* Full symmetry handling for linear molecules (otherwise up to D2h) 
* Parallelization using MPI library calls (MPI should be preinstalled) 
