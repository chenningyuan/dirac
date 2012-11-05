

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
