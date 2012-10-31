

.. _dirac:
.. index:: **DIRAC

==============================
\*\*DIRAC -- Job specification
==============================


Activation of modules
=====================

By default DIRAC only activates the generation of the basis set (e.g.  kinetic
balance conditions for small component functions) and the one-electron modules.
Two-electron integrals over the atomic basis functions will be calculated when
needed by the job modules given below.

We recommend that you use the `Dalton <http://www.daltonprogram.org>`_ program
package if you want to e.g. save the two-electron integrals to disk for another
purpose, this is not possible with DIRAC.


.. index:: .WAVE FUNCTION

.WAVE FUNCTION
--------------

Activates the wave function module(s).  This activates the reading of the
:ref:`wave_function` section, where the desired wave function type(s) must be
specified.  If you read in converged MO coefficients from the file DFCOEF and
you want to skip the SCF step and proceed directly to properties or the
post-SCF step, then a trick to do this is to comment out (or remove) this
keyword.


.. index:: .ANALYZE

.ANALYZE
--------

Activates the Hartree--Fock or Kohn--Sham analysis module.  This activates the
reading of the :ref:`analyze` section.


.. index:: .PROPERTIES

.PROPERTIES
-----------

Activates the property module (which will call the integral module for property
integrals).  This activates the reading of the :ref:`properties` section.


.. index:: .OPTIMIZE

.OPTIMIZE
---------

Activates the geometry optimization.  This activates the reading of the
:ref:`optimize` subsection.


.. index:: .4INDEX

.4INDEX
-------

Explicitly activates the transformation of integrals to molecular orbital
basis.  This activates the reading of the :ref:`moltra` section.

These transformed integrals are currently only used by the :ref:`relccsd`,
:ref:`reladc`, :ref:`dirrci`, and :ref:`lucita` modules, and if one of these
three modules are requested under :ref:`wave_function`, then this flag is
automatically activated unless .NO4INDEX is specified in this input module.

By default, molecular orbitals with orbital energy between -10 and +20 hartree
are included, this can be modified in the :ref:`moltra` section.


.. index:: .NO4INDEX

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


.. index:: .TITLE

.TITLE
------

Title line (max. 50 characters).
Example::

  .TITLE
   my first DIRAC calculation


.. index:: .INPTEST

.INPTEST
--------

Input test - no job modules are called, only verification of DIRAC input files.
It is often useful to start a new set of calculations with an input test in
order to check that input file processing is correct before submitting your
(long-term run) job.


.. index:: .ONLY INTEGRALS

.ONLY INTEGRALS
---------------

Stop after the calculation of the one-electron integrals for the Hamiltonian
and the one-electron integrals specified under :ref:`integrals`.  The integrals
are written to disk.


.. _optimize:
.. index:: *OPTIMIZE

\*OPTIMIZE
==========

Geometry optimization directives.

This section controls the geometry optimization. The geometry optimization
algorithm is based on the one of `Dalton <http://www.daltonprogram.org>`_. The
directives are therefore similar except for these few changes:

-  No second-order algorithms are available since the molecular Hessian
   is not implemented.
-  A few new keywords have been introduced (see below).

If :ref:`properties` or :ref:`analyze` is specified in the job input section
together with :ref:`optimize`, then the property or analysis module is called
in each optimization iteration and at the converged geometry.

Depending on convergence criteria for energies and gradients the values of the
thresholds may be automatically adjusted.


.. index:: .NO SKIP

.NO SKIP
--------

This keyword forces the LS and SS two-electron gradient to always be evaluated
in all geometry iterations (depending on the integral flag).

The SL and SS two-electron integral contributions to the gradient are normally
skipped if their contribution is estimated to be small. An "empirical" estimate
of the norms of the LS and SS two-electron gradients based on the norm of the
LL two-electron gradient.  This trick is by default activated in a geometry
optimization, since when the current geometry is far away from the equilibrium,
e.g. the norm of the gradient is, say, 1.0, then there is no need to calculate
the LS and/or SS two-electron gradient because they have a norm of, say, 0.001.


.. index:: .NUMGRA

.NUMGRA
-------

Force the use of numerical gradient in geometry optimization.


.. index:: .TWOGRD

.TWOGRD
-------

Include LL, SL, and SS integral contributions to the gradient
(1 = on; 0 = off). The default is to turn all on::

  .TWOGRD
   1 1 1


.. index:: .BAKER

.BAKER
------

Baker's convergence criteria [J. Comp. Chem. 14 (1993) 1085] will be
used.


.. index:: .1STORD

.1STORD
-------

Default 1st order method will be used: BFGS update.


.. index:: .GRADIENT

.GRADIENT
---------

Convergence threshold for the gradient. The default is 1.0D-5.


.. index:: .ENERGY

.ENERGY
-------

Convergence threshold for the energy. The default is 1.0D-6.


.. index:: .STEP THRESHOLD

.STEP THRESHOLD
---------------

Convergence threshold for the step. The default is 1.0D-5.
