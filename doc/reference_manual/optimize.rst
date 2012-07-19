

.. _optimize:

\*OPTIMIZE 
==========

Geometry optimization directives.

This section controls the geometry optimization. The geometry
optimization algorithm is based on the one of
`DALTON <http://www.kjemi.uio.no/software/dalton/dalton.html>`_. The
directives are therefore the same except for these few changes:

-  No second-order algorithms are available since the molecular Hessian
   is not implemented.
-  A few new keywords have been introduced (see below).

If :ref:`properties` or :ref:`analyze` is specified in the job input section
together with :ref:`optimize`, then the property or analysis module is called
in each optimization iteration and at the converged geometry.

Depending on convergence criteria for energies and gradients the values of the
thresholds may be automatically adjusted.


.NO SKIP
--------

The SL and SS two-electron integral contributions to the gradient are normally
skipped if their contribution is estimated to be small. An "empirical" estimate
of the norms of the LS and SS two-electron gradients based on the norm of the
LL two-electron gradient.  This trick is by default activated in a geometry
optimization, since when the current geometry is far away from the equilibrium,
e.g. the norm of the gradient is, say, 1.0, then there is no need to calculate
the LS and/or SS two-electron gradient because they have a norm of, say, 0.001.

This keywords forces the LS and SS two-electron gradient to always be evaluated
in all geometry iterations (depending on the integral flag).


.NUMGRA
-------

Force the use of numerical gradient in geometry optimization.


.TWOGRD
-------

Include LL, SL, and SS integral contributions to the gradient
(1 = on; 0 = off). Default::

  .TWOGRD
   1 1 1


.BAKER
------

Baker's convergence criteria [J. Comp. Chem. 14 (1993) 1085] will be
used.


.1STORD
-------

Default 1st order method will be used: BFGS update.


.GRADIENT
---------

Convergence threshold for the gradient. Default 1.0D-5.


.ENERGY
-------

Convergence threshold for the energy. Default 1.0D-6.


.STEP THRESHOLD
---------------

Convergence threshold for the step. Default 0.1D-4.
