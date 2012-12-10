

starstar(WAVE FUNCTION)


**Get the wave function**
=========================

This section allows the specification of which wave function module(s)
to activate. By default no modules are activated. To activate any of
these modules you must also specify :ref:`.WAVE FUNCTION`
under :ref:`**DIRAC`,
otherwise this input is not read.

Note that the order below specifies the order in which the different
modules are called if you ask for more than one.


keyword(SCF)

Activates the Hartree-Fock/Kohn-Sham module.

Specification of the SCF module can be given in the :ref:`*SCF` subsection.

If :ref:`.DFT` has been specified
under :ref:`**HAMILTONIAN`, then a
Kohn-Sham calculation will be performed, otherwise a Hartree-Fock
calculation will be performed.


keyword(RESOLVE)

Resolve open-shell states: do a small CI calculation to get the
individual energies of the states present in an
average-of-configurations open-shell Hartree-Fock calculation (see
:ref:`*RESOLVE`).


keyword(MP2)

Activates :ref:`.MP2CAL`.


keyword(MVO)

Calculate modified virtual orbitals (see :ref:`MVOCAL`).


keyword(RELCCSD)

Activates the  :ref:`RELCCSD`
(and the :ref:`MOLTRA`
module to get 4-index
transformed integrals).

By default, molecular orbitals with orbital energy between -10 and +20
hartree (a.u.) are included, this can be modified in the
:ref:`MOLTRA` section.


keyword(DIRRCI)

Activates the MOLFDIR CI module (and also the
:ref:`MOLTRA`  module to get 4-index transformed integrals).

Specification of input for the MOLFDIR CI module is given in
the :ref:`DIRRCI`  and :ref:`GOSCIP` sections.

By default, molecular orbitals with orbital energy between -10 and +20
hartree (a.u.) are included, this can be modified in the
:ref:`MOLTRA`  section.


keyword(LUCITA)

Activates the :ref:`*LUCITA` (and
the :ref:`**MOLTRA` module to get 4-index
transformed integrals).

By default, molecular orbitals with orbital energy between -10 and +20
hartree (a.u.) are included, this can be modified in the
:ref:`**MOLTRA` section.


**Pre-SCF orbital manipulations**
=================================

keyword(REORDER MO)

Interchange initial molecular orbitals prior to the SCF-calculation. The
start orbitals from DFCOEF are read and reordered.

For each fermion irrep give the new order of orbitals.

*Example:*

::

    .REORDER MO'S
    1..8,10,9


keyword(ORBROT)

Jacobi rotations between pairs of orbitals.

On the line following the keyword, give first the rotation angle, then
on the following line(s) for each fermion irrep, give an :ref:`orbital_strings`
of orbitals to
rotate.


**Post-SCF orbital manipulations**
==================================

keyword(POST SCF REORDER MO)

Interchange converged molecular orbitals. The orbitals from DFCOEF are
read and reordered just before exiting the SCF subroutine.

For each fermion irrep give the new order of orbitals.

*Example:*

::

    .POST DHF REORDER MO'S
    1..8,10,9


keyword(PHCOEF)

Phase adjustment of coefficients DFCOEF: make the largest element of a
given orbital real and positive.


keyword(RELADC)

Activates the :ref:`RELADC` and
calculates the ionization spectrum by the (A)lgebraic (D)iagrammatic
(C)onstruction ADC. Also activates the :ref:`**MOLTRA`
module to get 4-index
transformed integrals.


keyword(KRCI)

Activates the :ref:`KRCI` module
for the calculation of ground and excited states at the relativistic CI level.

keyword(KRMCSCF)

Activates the :ref:`KRMCSCF`
module for the optimization of ground and excited states (in other than the ground state symmetries) 
at the relativistic MCSCF level.


SCF & DFT
=========

This section gives directives for Hartree-Fock calculations.
Kohn-Sham calculations are performed by invoking the keyword
:ref:`.DFT` under :ref:`**HAMILTONIAN`.

Open-shell calculations correspond to either
average-of-configurations (:ref:`.AOC`) or fractional
occupation (:ref:`.FOCC`).
The former is the default for
Hartree-Fock calculations, whereas the latter is default for
Kohn-Sham calculations. Note that average-of-configurations
Kohn-Sham calculations are not well defined.

Occupation
==========

.. _.CLOSED SHELL:

.CLOSED SHELL
-------------

For each fermion irrep give the number of closed shell electrons.

The specification of the closed shell electrons is simple. For
symmetry groups *without* inversion symmetry, there is only one
fermion irrep, and you need only to specify the number of
electrons.

For symmetry groups *with* inversion symmetry, you need to specify
the distribution of the electrons in the two fermion irreps [Saue2000]_.

.. _.OPEN SHELL:

.OPEN SHELL
-----------

Specification of open shell(s).

For each open shell give the number of electrons and the number of
active spinors.

*Short example:*

::

    .OPEN SHELL
     1
     5/0,6

1 open shell with 5 electrons in 6 spinors (= 3 Kramers pairs) in
irrep 2. Thus, the fractional occupation is 5/6.

The open shell module in DIRAC is based on
average-of-configurations [Thyssen1998]_ .
The simplest case is one electron in two spinors (= one Kramers
pair). For this special case the average-of-configuration
calculation gives the same result as the usual restricted
open-shell Hartree-Fock. For all other cases the calculation gives
the average energy of many states.

Note that the order of closed and open shells are assumed to be as
in the following scheme:

**Virtual orbitals**
Not occupied

...
**Open shell 2**
Fractionally occupied

**Open shell 1**
Fractionally occupied

**Closed shell**

Doubly occupied;
that is, the lowest molecular orbitals are doubly occupied, the
next ones are occupied with the electrons of open shell no. 1,
etc.

Other orderings can be achieved by using :ref:`.REORDER MO'S`
and :ref:`.OVLSEL`.

To get the energies of the individual states present in the
average-of-configurations, specify :ref:`.RESOLVE` (see also
:ref:`*RESOLVE`).

To get the energies of (some) of the individual states present in
the average of configurations, you can use the
:ref:`GOSCIP`, the :ref:`DIRRCI`, or the :ref:`LUCITA`.

.BOSONS
-------

Occupation of boson irreps in spin-free calculation. For example,
for the D2h symmetry eight numbers in subsequent line, for the C2v
symmetry there four occupation numbers in line.

.MJSELE
-------

Occupation of MJ-splitted fermion irreps (only for the linear
symmetry). Next line contains sequence of occupation numbers, see
the SCF output for how many.

.. _.AOC:

.AOC
----

Average-of-configuration calculation (default for open-shell
Hartree-Fock).

.. _.FOCC:

.FOCC
-----

Fractional occupation (default for open-shell
Kohn-Sham) CLARIFY

:ref:`.FOCC` calculations are less memory-intensive than
:ref:`.AOC` calculations. In the latter case one additional
AO-Fock matrix is generated for each open shell.

:ref:`.FOCC` calculations are therefore an interesting option
for generating start orbitals for MCSCF as well as initial
convergence in open-shell Hartree-Fock.

.AUTOCC
-------

Program is allowed to change occupation during SCF cycles. This is
deactivated by default. However, the program will still try to do
an automatic initial occupation if neither
:ref:`.CLOSED SHELL` nor :ref:`.OPEN SHELL` is given.

Trial function
==============

An SCF-calculation (HF or DFT) may be initiated in three different
ways:


-  Using **MO coefficients** from a previous calculation.
-  Using coefficients obtained by diagonalization of the
   one-electron Fock matrix: the **bare nucleus approach**.
-  Using the **two-electron Fock matrix** from a previous
   calculation; this may be thought of as starting from a converged
   SCF potential.

The default is to start from MO coefficients if the file DFCOEF is
present. Otherwise the corrected bare nucleus approach is followed.
In all three cases linear dependencies are removed in the zeroth
iteration.

.NOBNCR
-------

Switch off the bare nucleus correction. This correction is on per
default and improves the screening of the nuclei by estimating
two-electron repulsion via nuclear-attraction type integrals:

.. math::

 \langle X_{A} \vert \sum_{C} \frac{-Z_{C} \cdot \sum_j a_j e^{(-\alpha^C_j r_{C}^{2})}}{r_{C}} \vert X_{B} \rangle, \ \ \ X = L,S

The coefficients *a*\ :sub:`*j*`\  and the exponents :math:`\alpha^{C_j}`
in this expression are chosen according to Slater's rules to obtain
an approximate atomic electronic density for the initial guess. For
example, with one heavy element and without this correction (that
is, with the bare nucleus Hamiltonian) all electrons will end up on
that heavy element in the initial guess!

.TRIVEC
-------

Start SCF-iterations from the vector file.

.TRIFCK
-------

Start SCF-iterations from the two-electron Fock matrix from
previous calculation (stored on file DFFCK2).

Convergence criteria
====================

Three different criteria for convergence may be chosen:


-  The norm of the DIIS error vector
   :math:`\mathbf{e} = [\mathbf{F}, \mathbf{D}]` (in MO basis). This
   corresponds to the norm of the electronic gradient and is the
   recommended convergence criterion. When you are only interested in
   the energy :ref:`.EVCCNV` = 1.0D-5 is usually sufficient.
   For properties and correlated methods you should converge to
   :ref:`.EVCCNV` = 1.0D-9. Large negative energy eigenvalues
   lead to a loss of precision that might lead to convergence
   problems. Remember also that a too loose screening threshold (too
   many integrals neglected) will hinder convergence. You should
   modify
   :ref:`.SCREEN` under
   :ref:`*TWOINT` if you
   modify :ref:`.EVCCNV` or one of the other two convergence
   criteria.
-  The difference in total energy between two consecutive
   iterations.
-  The largest absolute difference in the total Fock matrix between
   two consecutive iterations.

The change in total energy is approximately the square of the
largest element in the error vector or the largest change in the
Fock matrix. The default is convergence on electronic gradient with
1.0D-6 as threshold. Alternatively, the iterations will stop at the
maximum number of iterations.

Sometimes it may happen that the specified convergence criterion is
too tight for the given basis set and/or other input parameters. In
this case one needs to decide whether one should proceed with
post-HF steps (like correlation calculations) or not. The program
decides this by looking at a secondary convergence criterion that
gives the *allowed* convergence. This value is by default the same
as first or *desired* convergence criterion but can be made lower
to make sure that a calculation does not abort when the convergence
is slightly above threshold.

For more detailed help see SCF help on convergence troubleshooting  and related pages.

.MAXITR
-------

Maximum number of SCF iterations.

*Default:*

::

    .MAXITR
     50

When restarting SCF itrations from previous molecular orbitals file
(DFCMO or formatted DFPCMO), we recommend to decrease maximum
number of iterations together with readjusting desired and allowed
convergence thresholds. By properly set desired and allowed
thresholds one can have exact number of iterations specified by
*.MAXITR*.

.. _.EVCCNV:

.EVCCNV
-------

Converge on error vector (electronic gradient).

*2 (real) Arguments:*

::

    .EVCCNV
     desired threshold allowed threshold

.. _.ERGCNV:

.ERGCNV
-------

Threshold for convergence on total energy.

*2 (real) Arguments:*

::

    .ERGCNV
     desired threshold allowed threshold

.. _.FCKCNV:

.FCKCNV
-------

Converge on largest absolute change in Fock matrix.

*2 (real) Arguments:*

::

    .FCKCNV
     desired threshold allowed threshold

Note that the *allowed* threshold may be omitted. It is then made
equal to the *desired* threshold.

Convergence acceleration
========================

It is imperative to keep the number of SCF iterations at a minimum.
This may be achieved by convergence acceleration schemes:


-  **Damping:** The simplest scheme is damping of the Fock matrix
   that may remove oscillations. In :math:`n + 1` iteration the Fock matrix
   to be diagonalized is:
   :math:`\mathbf{F}\' = (1-c) \mathbf{F}_{n+1} + c \mathbf{F}_n`,
   where :math:`c` is the damping factor.

-  **DIIS:** Direct inversion of iterative subspaces, Refs. [Pulay1980]_ , [Pulay1982]_ , [Hamilton1986]_,
   may be thought of as generalized damping involving Fock matrices from many iterations.
   Damping factors are obtained by solving a simple matrix equation
   involving the B-matrix constructed from error vectors (approximate
   gradients). Linear dependent columns in the B-matrix is removed.

In DIRAC DIIS takes precedence over damping.

.DIISTH
-------

Change the default convergence threshold for initiation of DIIS,
based on largest element of error vector.

*Default:*

::

    .DIISTH
     a very large number

.MXDIIS
-------

Maximum dimension of B-matrix in the DIIS module.

*Default:*

::

    .MXDIIS
     10

.. _.DIISMO:

.DIISMO
-------

Activate DIIS in orthogonal basis (MO) with the error vector as
described above.

.. _.DIISAO:

.DIISAO
-------

Activate
`DALTON <http://www.kjemi.uio.no/software/dalton/dalton.html>`_-like
DIIS using AO-basis. The error vector is

.. math::

 \mathbf{e} = \mathbf{f} - \mathbf{f}^{\dagger}|

where the :math:`\mathbf{f}` is given by

.. math::

   \mathbf{f}=\mathbf{C}^{\dagger}\cdot\mathbf{S}_{AO}\cdot \left[ \mathbf{D}^{C}_{AO}\cdot\mathbf{F}^{D}_{AO}+
   \sum_{O\in\mathcal{O}}f_{O}\cdot\mathbf{D}^{O}_{{AO}}\cdot ( \mathbf{F}^{D}_{{AO}}+
   (a_{O}-1)\mathbf{Q}^{V,O}_{{AO}} ) \right] \mathbf{C}

.. _.NODIIS:

.NODIIS
-------

Do not perform DIIS. The default is to activate
:ref:`.DIISMO` for closed-shell calculations, and to
activate :ref:`.DIISAO` for average-of-configurations
calculations.

.DAMPFC
-------

Change the default damping factor.

*Default:*

::

    .DAMPFC
     0.25

.NODAMP
-------

Do not perform damping of the Fock matrix. Damping is activated by
default, but DIIS takes precedence. In case all columns in the
B-matrix is removed by linear dependency, damping is activated.

Level shift
===========

.LSHIFT
-------

Activate level shift (for virtual orbitals). Followed by a real
argument (level shift).

.OLEVEL
-------

Activate level shift (for open-shell orbitals). Followed by a real
argument (level shift), one line for each open shell {Please give
example}.

State selection
===============

Convergence can be improved by selection of vectors based on
overlap with vectors from a previous iteration. This method may
also be used for convergence to some excited state.

If dynamic overlap selection is used, the vector set from the
previous iteration is used as the criterion. For the first
iteration either restart vectors or vectors generated by the bare
nucleus approach (not*recommended) are used.*

If :ref:`.NODYNSEL` is given, either the restart vectors
or the bare nucleus vectors are used, i.e. the overlap selection
vectors are *not* updated in each iteration. Please note, that
overlap selection based on vectors from the bare nucleus approach
is not recommended.

Overlap selection is very useful together with
:ref:`.REORDER MO'S`.
This will reshuffle the vectors within the restart coefficients.

Example: First one might do a open shell calculation on Boron, this
would give the *P*\ :sub:`1 / 2`\  state. But if we restart on the
*P*\ :sub:`1 / 2`\  coefficients, interchange the
*p*\ :sub:`1 / 2`\  with the *p*\ :sub:`3 / 2`\  orbitals, and
request overlap selection, we can converge to the
*P*\ :sub:`3 / 2`\  state.

There also exists a keyword for reordering the converged SCF
orbitals. This is useful for reordering the orbitals for the
4-index transformation and subsequent correlation calculations
(CCSD, CI etc.) (see :ref:`.POST DHF REORDER MO'S`).

.. _.OVLSEL:

.OVLSEL
-------

Activate dynamic overlap selection. The default is no overlap
selection.

.. _.NODYNSEL:

.NODYNSEL
---------

No dynamic update of overlap selection vectors. The default is
dynamic update.

Iteration speedup
=================

The total run time may be reduced significantly by reducing the
number of integrals to be processed in each iteration:


-  **Screening on integrals:** Thresholds may be set to eliminate
   integrals below the threshold value, see [Saue1997]_.
   . The threshold for LL
   integrals is set in the basis file, but this threshold may be
   adjusted for SL and SS integrals by threshold factors set in the
   :ref:`**INTEGRALS` section.

-  **Screening on density:** In direct mode further reductions are
   obtained by screening on the density matrix as well, see Ref. [Saue1997]_.
   This becomes even more
   effective if one employs *differential densities*, that is
   :math:`\Delta \mathbf{D} = \mathbf{D}_{n+1} - \alpha \cdot \mathbf{D}_n`.
   The default value for :math:`alpha` is
   :math:`\alpha=\frac{ \mathbf{D}_{n+1} \cdot \mathbf{D}_n }{ \mathbf{D}_{n}  \cdot \mathbf{D}_n }`
   which corresponds to a Gram-Schmidt orthogonalization. As SCF
   converges, :math:`alpha` goes towards 1, but  :math:`alpha` can also explicitly be set equal
   to 1 with :ref:`.FIXDIF`.


-  **Neglect of integrals:** The number of integrals to be
   processed may be reduced even further by adding SL and SS integrals
   only at an advanced stage in the SCF-iterations, as determined
   either by the number of iterations or by energy convergence. The
   latter takes precedence over the former.

.NODSCF
-------

Do not perform SCF-iterations with differential density matrix.

*Default:* Use differential density matrix in direct SCF.

.. _.FIXDIF:

.FIXDIF
-------

Set :math:`alpha` equal to 1.

.CNVINT
-------

Set threshold for convergence before adding SL and SS integrals to
SCF-iterations.

*2 (real) Arguments:*

::

    .CNVINT
     CNVINT(1) CNVINT(2)

.ITRINT
-------

Set the number of iterations before adding SL and SS integrals to
SCF-iterations.

*Default:*

::

    .ITRINT
     1 1

.INTFLG
-------

Specify what two-electron integrals to include (see
:ref:`.INTFLG` under :ref:`**HAMILTONIAN`).

*Default:* :ref:`.INTFLG` from :ref:`**HAMILTONIAN`.

Output control
==============

.PRINT
------

General print level for the SCF method. For instance, value of 2
prints eignevalues during each iteration.

*Default:*

::

    .PRINT
     0

.EIGPRI
-------

Controls the print-out of positive energy and negative energy
eigenvalues (1 = on; 0 = off).

*Default:* Only the positive energy eigenvalues are printed.

::

    .EIGPRI
     1 0

Eliminating/freezing orbitals
=============================

In studies of electronic structure it may be of interest to
eliminate or freeze certain orbitals. This option is furthermore
useful for convergence, in particular to excited electronic states.
A simple case is the thallium atom. The ground state
\ :sup:`2`\ *P*\ :sub:`1 / 2`\  has the electronic configuration
[Xe]|4f^{14}5d^{10}6s^{2}6p^1\_{1/2}|. The first excited state
\ :sup:`2`\ *P*\ :sub:`3 / 2`\  with the electronic configuration
[Xe]4f^{14}5d^{10}6s^{2}6p^1\_{3/2}| can easily be accessed by
first calculating the ground state, then eliminating the
6p^1\_{1/2}| from the ensuing calculation. In the final
calculation the excited state is relaxed using overlap selection.
The use of frozen orbitals is demonstrated in test/33.frozen: When
the geometry of the water molecule is optimized with the oxygen 1s
and 2s orbitals frozen, a bond angle of 96.242 degrees is found,
contrary to the 90 degrees one might have expected when s-p
hybridization is thus blocked [Kutzelnigg1984]_.

The elimination of orbitals is achieved by projecting the selected
orbitals out of the transformation matrix to orthonormal basis. The
selected orbitals can be expressed either in the full molecular
basis or in the basis set of the chosen fragment. In the latter
case, the same set of fragment orbitals can in the case of atomic
fragments be used at different molecular geometries. One may even
perform a geometry optimization, but only using the numerical
gradient. When freezing orbitals the selected orbitals are first
eliminated from the transformation to orthonormal basis, but then
reintroduced in the backtransformation step. They will appear in
the output with zero orbital eigenvalues. Note that when freezing
orbitals the orbitals to be eliminated must be specified as well.
The frozen orbitals must be a subset of the eliminated orbitals.

Fragments are defined with respect to the list of
symmetry-independent atoms appearing in the DIRAC basis file:
Consider the water molecule in the full *C*\ :sub:`2*v*`\
symmetry. Then there are two fragments: the oxygen atom and the
H\ :sub:`2`\  moiety. However, with no symmetry there will be three
fragments: the oxygen atom and the two hydrogen atoms. At the
moment there are no orthonormalization of fragments on different
fragments and so in practice one should only use orbitals from one
fragment.

.PROJECT
--------

Eliminate orbitals by projecting them out of the transformation
matrix to orthonormal basis.

*Arguments:* Number of fragments (NPRJREF).

Then for each fragment read name PRJFIL of the coefficient file
followed by the number of symmetry-independent nuclei in this
fragment followed by an
:ref:`orbital_strings` of
selected orbitals for each fermion irrep.

::

    DO J = 1,NPRJREF
      READ(LUCMD,'(A6)') PRJFIL(J)
      READ(LUCMD,'(I6)') NPRJNUC(J)
      READ(LUCMD,'(A72)') (VCPROJ(I,J),I=1,NFSYM)
    ENDDO

.PRJTHR
-------

Smallest norm accepted when eliminating orbitals.

*Default:*

::

    .PRJTHR
     1.0D-10

.OWNBAS
-------

Eliminated/frozen orbitals are given in the fragment basis. Note
that the list of fragments is assumed to follow the list of
symmetry independent nuclei in the DIRAC basis file.

.FROZEN
-------

Freeze orbitals (must only be used in conjunction with
:ref:`.PROJECT`).

For each fermion irrep, give an :ref:`orbital_strings` of
orbitals to freeze.

.. _resolve:

\*RESOLVE
=========

Resolve open-shell states.

The individual electronic states of an average-of-configurations
calculation are resolved by a full CI in the small open-shell
manifold(s) using the :ref:`GOSCIP` module.

Note that even if this is a small CI diagonalization one needs in
principle to generate *all* AO integrals for the initial 4-index
transformation. The cost of calculation can be greatly reduced invoking
integral screening (controlled by :ref:`.SCREEN`) in the 4-index
transformation, possibly also by leaving out at least SS integrals
(using :ref:`.INTFLG`).

.PRINT
------

Print level.

*Default:*

::

    .PRINT
     0

.INTFLG
-------

Specify what two-electron integrals to include (default: :ref:`.INTFLG` under :ref:`**HAMILTONIAN`).

.SCREEN
-------

Screening threshold in the 4-index transformation.

A negative number deactivates screening.

*Default:*

::

    .SCREEN
     1.0D-14

.SCHEME
-------

Choose integral transformation algorithm (PLEASE CLARIFY!)

*Default:*

::

    .SCHEME
     4

\*MP2CAL -- MP2 calculation
===========================

This section gives directives for the integral-direct closed-shell MP2
calculation. The default algorithm is described
`here <http://folk.uio.no/vebjornb/usit/dirac-mp2/dirac-mp2.php>`_ and
is a modification of the original algorithm published in Ref. [Laerdahl1997]_.

.OCCUPIED
---------

Active, doubly occupied, electronic solutions.

For each fermion irrep, give an :ref:`orbital_strings`
of active orbitals.

*Default:* The occupied electronic solutions.

.VIRTUAL
--------

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


star(MVOCAL)

Modified virtual orbitals.

The virtual canonical Hartree-Fock orbitals are not the optimal choice
for correlated calculations in which the virtual space is truncated; the
lower orbitals tend to be diffuse since they see the *N* + 1 electron
system. One may exploit the freedom of rotation amongst the virtual
orbitals to obtain orbitals better suited for correlated calculations.
The default option in DIRAC is to generate the modified virtual orbitals
by construction of a Fock operator for the system with some electrons
removed, as first proposed by Bauschlicher [Bauschlicher1980]_.
The orbital string
specified with :ref:`.MVOVEC` indicates what occupied orbitals to
remove in the construction of the Hartree-Fock potential.

*Note* that the resulting modified virtual orbitals are no longer
canonical and can *not* be used in MP2 calculations; they can be used in
the :ref:`DIRRCI` and
:ref:`RELCCSD` modules.


keyword(PRINT)

Print level.

*Default:*

::

    .PRINT
     0


keyword(VECMVO)

Read :ref:`orbital_strings` of
orbitals to be included in the Hartree-Fock potential when constructing
the cationic Fock operator.


keyword(INTFLG)

Specify what two-electron integrals to include during the construction
of the modified virtual orbitals (default: :ref:`.INTFLG` under :ref:`**HAMILTONIAN`).
