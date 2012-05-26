.. _dirac:

\*\*DIRAC -- Job specification
==============================

By default DIRAC only activates the generation of the basis set (e.g.
kinetic balance conditions for small component functions) and the
one-electron modules. Two-electron integrals over the atomic basis
functions will be calculated when needed by the job modules given below.

We recommend that you use the
`DALTON <http://www.kjemi.uio.no/software/dalton/dalton.html>`_ program
package if you want to e.g. save the two-electron integrals to disk for
another purpose, this is not possible with DIRAC.


`.WAVE FUNCTION`_
-----------------

Activates the wave function module(s).  This activates the reading of the
:ref:`wave_function` section, where the desired wave function type(s) must be
specified.


.ANALYZE
--------

Activates the Hartree-Fock or Kohn-Sham analysis module.  This activates the
reading of the :ref:`analyze` section.


.PROPERTIES
-----------

Activates the property module (which will call the integral module for property
integrals).  This activates the reading of the :ref:`properties` section.


.OPTIMIZE
---------

Activates the geometry optimization.  This activates the reading of the
:ref:`optimize` subsection.


.4INDEX
-------

Explicitly activates transformation of integrals to molecular orbital basis.
This activates the reading of the :ref:`moltra` section.

These transformed integrals are currently only used by the
:ref:`relccsd`,
:ref:`reladc`,
:ref:`dirrci`, and
:ref:`lucita`
modules,
and if one of these three modules
are requested under :ref:`wave_function`, then this flag is automatically
activated unless .NO4INDEX is specified in this input module.

By default, molecular orbitals with orbital energy between -10 and +20 hartree
(a.u.) are included, this can be modified in the :ref:`moltra` section.


.NO4INDEX
---------

Do not automatically activate integral transformation to molecular
orbital basis if any of
:ref:`relccsd`,
:ref:`reladc`,
:ref:`dirrci`, and
:ref:`lucita`
modules
are requested under :ref:`wave_function`.

This keyword is utilized when repeating correlated CC or CI calculations
(with different parameters for instance) based on saved files after the
integral transformation.


.TITLE
------

Title line (max. 50 characters). Example::

  .TITLE
   my first DIRAC calculation


.INPTEST
--------

Input test - no job modules are called, only verification of DIRAC input files.
It is often useful to start a new set of calculations with an input test in
order to check that input file processing is correct before submitting your
(long-term run) job.


.ONLY INTEGRALS
---------------

Stop after the calculation of the one-electron integrals for the Hamiltonian
and the one-electron integrals specified under :ref:`integrals`.  The integrals
are written to disk.
