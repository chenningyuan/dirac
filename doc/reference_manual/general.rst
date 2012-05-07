Manual:GENERAL
==============

\*\*GENERAL -- General input
============================

..`|image0| </dirac/index.php/Image:Dirac.png>`_

**Advanced options**

.DIRECT
-------

Direct evaluation of two-electron integrals for Fock matrices (all
two-electron integrals for other uses, e.g. CI, CCSD, MCSCF, are always
evaluated directly, i.e. never read from disk).

*Default:* Evaluate LL, SL, and SS integrals directly (1 = evaluate
directly; 0 = do not evaluate directly).

::

    .DIRECT
     1 1 1

.SPHTRA
-------

Transformation to spherical harmonics embedded in the transformation to
orthonormal basis; totally symmetric contributions are deleted.

*Default:*

::

    .SPHTRA
     1 1

Spherical transformation of large and small components, respectively (1
= on; 0 = off).

The transformation of the large components is a standard transformation
from Cartesian functions to spherical harmonics. The transformation of
the small components is modified, however, in accordance with the
restricted kinetic balance relation.

.PCMOIN
-------

Read MO coefficients from the formatted file DFPCMO.

This is useful for porting coefficients between machines with different
binary structure, for example, between Silicon Graphics and IBM. See
also `.PCMOUT <#.PCMOUT>`_.

.PCMOUT
-------

Write MO coefficients to the formatted file DFPCMO.

This is useful for porting coefficients between machines with different
binary structure, for example, between Silicon Graphics and IBM. See
also `.PCMOIN <#.PCMOIN>`_.

.ACMOUT
-------

Dump coefficients in *C*\ :sub:`1`\  symmetry to the unformatted file
DFACMO.

This can be useful, e.g. when doing a projection analysis: The molecule
is calculated at its highest symmetry, each constituent atom type is
calculated in linear symmetry (specifying no symmetry operations in the
`basis file </dirac/index.php/Manual:Basis:Title_and_symmetry>`_), and
then the analysis is carried out in *C*\ :sub:`1`\  symmetry, or the
coefficients are adapted to current symmetry (see
`.ACMOIN <#.ACMOIN>`_).

.ACMOIN
-------

Import coefficients in *C*\ :sub:`1`\  symmetry from unformatted file
DFACMO to current symmetry. This assumes that the current symmetry is
lower than the symmetry used for obtaining the original coefficients.

.LOWJACO
--------

Use Jacobi diagonalization in the Löwdin orthogonalization procedure.
This is much slower than the default method but does not mix AOs in the
case of block-diagonal overlap matrix. This keyword is only available in
the development version, DIRAC 08 uses Jacobi diagonalization by
default.

.DOJACO
-------

The default method is to use the
`Householder <http://www.maths.lancs.ac.uk/~gilbert/m306c/node21.html>`_
diagonalization method, which is generally more efficient. This method
may, however, mix degenerate eigenvectors of different symmetries. The
Jacobi method is currently limited to real groups.

.QJACO
------

Employ pure Jacobi diagonalization of quaternion matrixes. Properly
handles degenerate eigenvectors. Slower than .DOJACO and exclusive.
Experimental option.

.LINDEP
-------

Thresholds for linear dependence in large and small components; refer to
the smallest acceptable values of eigenvalues of the overlap matrix.

*Default:*

::

    .LINDEP
     1.0D-6 1.0D-8

.RKBIMP
-------

Read in SCF coefficients calculated using restricted kinetic balance
(RKB) and add the UKB component (all small component). This option is a
nice way to get (unrestricted) magnetic balance in response calcualtions
involving a uniform magnetic field (e.g. NMR shielding and
magnetizability), in particular when combined with London orbitals,
which makes the magnetic balance *atomic*.

.PRJTHR
-------

RKBIMP projects out the RKB coefficients transformed to orthonormal
basis and then adds the remained, corresponding to the UKB complement.
With the keyword you can set the threshold for projection.

*Default:*

::

    .PRJTHR
    1.0D-5

`|image1| </dirac/index.php/Image:HJAaJ.png>`_

**Programmers options**

.PRINT
------

General print level. The higher the number is, the more output the user
gets. Option of this type is useful for code debugging.

*Default:*

::

    .PRINT
     0

.CVALUE
-------

Speed of light in a.u.

*Default:*

::

    .CVALUE
     137.0359998

.LOGMEM
-------

Write out a line in the output for each memory allocation done by Dirac.
This is mainly useful for programmers or for solving out-of-memory
issues.


\*PARALLEL -- Parallelization directives
========================================

.. `|image0| </dirac/index.php/Image:HJAaJ.png>`_

**Programmers options**

.PRINT
------

Print level for the parallel calculation.

*Default:*

::

    .PRINT
     0

A print level of at least 2 is needed in order to be able to evaluate
the parallelization efficiency. A complete timing for all nodes will be
given if the print level is 4 or higher.

.NTASK
------

Number of tasks to send to each node when distributing the calculation
of two-electron integrals.

*Default:*

::

    .NTASK
     1

A task is defined as a shell of atomic integrals, a shell being an input
block. One may therefore increase the number of shells given to each
node in order to reduce the amount of communication. However, the
program uses dynamical allocation of work to each node, and thus this
option should be used with some care, as too large tasks may cause the
dynamical load balancing to fail, giving an overall decrease in
efficiency. The parallelization is also very coarse grained, so that the
amount of communication seldom represents any significant problem.

Sample inputs
=============

Use a reduced value for the speed of light to amplify the effect of
relativity (for testing):

::

    **GENERAL
    .CVALUE
     50.0

Very often one needs to export or to read formatted file with orbital
coefficients - DFPCMO. Therefore the user should provide for the reading

::

    **GENERAL
    .PCMOIN

or

::

    **GENERAL
    .PCMOUT

for the export.

