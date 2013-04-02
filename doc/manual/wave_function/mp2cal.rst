star(MP2CAL)

**Introduction**
================

This section gives directives for the integral-direct closed-shell 
second-order Møller-Plesset perturbation theory (MP2)
calculation. The default algorithm is described
`here <http://folk.uio.no/vebjornb/usit/dirac-mp2/dirac-mp2.php>`_ and
is a modification of the original algorithm published in Ref. [Laerdahl1997]_.

**Basic keywords**
==================

keyword(OCCUPIED)

Active, doubly occupied, electronic solutions.

For each fermion irrep, give an :ref:`orbital_strings`
of active orbitals.

*Default:* The occupied electronic solutions.

keyword(VIRTUAL)

Active virtual solutions.

For each fermion irrep, give an `orbital_strings`
of active orbitals.

*Default:* All the unoccupied virtual solutions.


**Advanced options**
====================

keyword(INTFLG)

Specify what two-electron integrals to include in the direct MP2
calculation (default: :ref:`.INTFLG` under :ref:`**HAMILTONIAN`).


keyword(SCREEN)

Screening threshold in the 4-index transformation.

A negative number deactivates screening.

*Default:*

::

    .SCREEN
     1.0D-14


**Programmers options**
=======================

keyword(PRINT)

Print level.

*Default:*

::

    .PRINT
     0


keyword(VIRTHR)

Threshold for neglecting high-lying virtual orbitals. Neglect all
orbitals with energy above this threshold. This option is retained for
backwards compatibility, it is advised to use the more general energy
option in the :ref:`orbital_strings`
input.

*Arguments:* Real DMP2\_VIRTHR.

*Default:* No threshold.


keyword(IJTSK)

Max. number of I (or J) orbitals (occupied orbitals) in each calculation
batch. Usually one batch is sufficient but larger calculations may be
possible by reducing the batch size. Reducing the batch size to half
will approximately reduce the memory requirements to 1/4 and increase
the CPU time by a factor 4.

*Arguments:* Integer IJTSK.

*Default:* All active occupied electronic solutions.


keyword(SCLMEM)

Internal memory available for scalar untransformed integrals. Larger
batches will be written to disk.

*Arguments:* Integer MAXSCL.

*Default:* Set by program and is usually sufficient.


keyword(ORGALG)

Use original alorithm as described in Ref. [Laerdahl1997]_.

*Default:* Use `new
algorithm <http://folk.uio.no/vebjornb/usit/dirac-mp2/dirac-mp2.php>`_.


