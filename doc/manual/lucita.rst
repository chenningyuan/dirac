

star(LUCITA)

Spin-free relativistic GASCI module written by Jeppe Olsen, adaptation to DIRAC by Timo Fleig, parallelization by Stefan
Knecht and Hans Joergen Aagaard Jensen.

LUCITA is a string-based Hamiltonian-direct configuration interaction (CI) program, based on the LUCIA code [Olsen1990]_.
It is capable of doing efficient CI computations at arbitrary excitation level, e.g. FCI, SDCI, RASCI, and MRCI using
general active spaces. The code is interfaced [Fleig2005]_ to molecular integrals obtained in the spin-free Dirac
formalism and uses non-relativistic point group symmetry. It is implemented as a full parallel version [Knecht2008]_.

A central feature of the program is the Generalized Active Space (GAS) concept, in which the underlying total orbital
space is subdivided into a basically arbitrary number (save for an upper limit) of subspaces with arbitrary occupation
constraints. This is the most general approach to orbital space subdivisions. The program uses DIRAC orbitals from either
a closed- or an open-shell calculation in a spin-free relativistic formalism or the non-relativistic Lévy-Leblond
formalism (see :ref:`**HAMILTONIAN`).

The technical limitations are roughly set by several 100 million determinants in the CI expansion on PCs and common computing clusters and several billions of determinants on supercomputers with ample memory.
The program also computes 1- and, if requested 2-particle densities from optimized CI wave functions. The density matrices may be printed along with natural orbital occupations and the corresponding eigenvectors (NOs).

**General information**
=======================

Please add the lines

::

    .SCHEME
     4

to the :ref:`**MOLTRA` input section for the integral transformation.

**Mandatory keywords**
======================

keyword(INIWFC)

specifies the initial SCF wave function. 

*closed-shell:*

::

    .INIWFC
     DHFSCF

*open-shell:*

::

    .INIWFC
     OSHSCF

keyword(CITYPE)

Type of CI calculation. Typical multi-reference CI calculations should be defined by RASCI or GASCI.

valid choices are: FCI, SDCI, SDTQ, RASCI, GASCI.

*example:*

::

    .CITYPE
     FCI

keyword(MULTIP)

specifies the target state spin multiplicity (2S+1) of the desired eigenstates.

*example for singlet state(s) (S=0):*

::

    .MULTIP
     1

