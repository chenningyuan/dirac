:orphan:
 

star(ESR)

This section gives directives for the calculation of ESR/EPR parameters, g-tensors and Hyperfine Coupling tensors (HFC).

..
  COMMENT. Following lines are from nmr.rst, however, LONDON is not implemented for *ESR yet. /Oct 2022
  If common gauge origin (CGO) is used, i.e. :LONDON not specified,
  then the user can define the gauge origin used for the external
  magnetic field with :ref:`HAMILTONIAN_.GAUGEORIGIN` or :ref:`HAMILTONIAN_.GO ANG`
  under :ref:`**HAMILTONIAN`.
  If :ref:`NMR_.USECM` is specified then center-of-mass is used as gauge origin.
  Default is (0, 0, 0).


General control statements
==========================

keyword(PRINT)

Print level.

*Default:*

::

    .PRINT
     0

Definition of the multiplicity and properties
=============================================

keyword(MULTIP)

*Default:*

Doublet for an odd number of electrons, triplet for an even number of electrons.

keyword(STATE)

The first CI root belonging to the multiplet of quasi-degenerate states.
For example, if the ground state is a non-degenerate singlet state (non-relativistically a spin singlet)
followed by the triplet state you want to study, then you should specify state no. 2 here.

*Default:*

::

    .STATE
     1

keyword(SELECT)

Select which nuclei to calculate HFC for. Default is *all*, which can be quite time-consumig if the molecule has many atoms.
(You cannot use :ref:`INTEGRALS_.SELECT` under :ref:`**INTEGRALS` if you use this keyword, they do the same thing.
It is included here for convenience.)
To select three nuclei, no. 3, 7, and 8:

::

     .SELECT
      3
      3 7 8

..
  **Advanced options**

