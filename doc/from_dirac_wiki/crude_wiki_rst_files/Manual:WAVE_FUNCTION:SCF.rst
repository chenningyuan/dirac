Manual:WAVE FUNCTION:SCF
========================

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

Contents
--------

-  `1 Occupation <#Occupation>`_

   -  `1.1 .CLOSED SHELL <#.CLOSED_SHELL>`_
   -  `1.2 .OPEN SHELL <#.OPEN_SHELL>`_
   -  `1.3 .BOSONS <#.BOSONS>`_
   -  `1.4 .MJSELE <#.MJSELE>`_
   -  `1.5 .AOC <#.AOC>`_
   -  `1.6 .FOCC <#.FOCC>`_
   -  `1.7 .AUTOCC <#.AUTOCC>`_

-  `2 Trial function <#Trial_function>`_

   -  `2.1 .NOBNCR <#.NOBNCR>`_
   -  `2.2 .TRIVEC <#.TRIVEC>`_
   -  `2.3 .TRIFCK <#.TRIFCK>`_

-  `3 Convergence criteria <#Convergence_criteria>`_

   -  `3.1 .MAXITR <#.MAXITR>`_
   -  `3.2 .EVCCNV <#.EVCCNV>`_
   -  `3.3 .ERGCNV <#.ERGCNV>`_
   -  `3.4 .FCKCNV <#.FCKCNV>`_

-  `4 Convergence acceleration <#Convergence_acceleration>`_

   -  `4.1 .DIISTH <#.DIISTH>`_
   -  `4.2 .MXDIIS <#.MXDIIS>`_
   -  `4.3 .DIISMO <#.DIISMO>`_
   -  `4.4 .DIISAO <#.DIISAO>`_
   -  `4.5 .NODIIS <#.NODIIS>`_
   -  `4.6 .DAMPFC <#.DAMPFC>`_
   -  `4.7 .NODAMP <#.NODAMP>`_

-  `5 Level shift <#Level_shift>`_

   -  `5.1 .LSHIFT <#.LSHIFT>`_
   -  `5.2 .OLEVEL <#.OLEVEL>`_

-  `6 State selection <#State_selection>`_

   -  `6.1 .OVLSEL <#.OVLSEL>`_
   -  `6.2 .NODYNSEL <#.NODYNSEL>`_

-  `7 Iteration speedup <#Iteration_speedup>`_

   -  `7.1 .NODSCF <#.NODSCF>`_
   -  `7.2 .FIXDIF <#.FIXDIF>`_
   -  `7.3 .CNVINT <#.CNVINT>`_
   -  `7.4 .ITRINT <#.ITRINT>`_
   -  `7.5 .INTFLG <#.INTFLG>`_

-  `8 Output control <#Output_control>`_

   -  `8.1 .PRINT <#.PRINT>`_
   -  `8.2 .EIGPRI <#.EIGPRI>`_

-  `9 Eliminating/freezing orbitals <#Eliminating.2Ffreezing_orbitals>`_

   -  `9.1 .PROJECT <#.PROJECT>`_
   -  `9.2 .PRJTHR <#.PRJTHR>`_
   -  `9.3 .OWNBAS <#.OWNBAS>`_
   -  `9.4 .FROZEN <#.FROZEN>`_

-  `10 Sample input <#Sample_input>`_
-  `11 References <#References>`_

This section gives directives for Hartree-Fock calculations. Kohn-Sham
calculations are performed by invoking the keyword
`.DFT </dirac/index.php/Manual:HAMILTONIAN#.DFT>`_ under
`\*\*HAMILTONIAN </dirac/index.php/Manual:HAMILTONIAN>`_. Open-shell
calculations correspond to either average-of-configurations
(`.AOC <#.AOC>`_) or fractional occupation (`.FOCC <#.FOCC>`_). The
former is the default for Hartree-Fock calculations, whereas the latter
is default for Kohn-Sham calculations. Note that
average-of-configurations Kohn-Sham calculations are not well defined.

Occupation
==========

.CLOSED SHELL
-------------

For each fermion irrep give the number of closed shell electrons.

The specification of the closed shell electrons is simple. For symmetry
groups *without* inversion symmetry, there is only one fermion irrep,
and you need only to specify the number of electrons.

*Example:* Minimum wave function input for water (the *C*\ :sub:`2*v*`\ 
point group):

::

    **WAVE FUNCTION
    .SCF
    *SCF
    .CLOSED SHELL
     10

For symmetry groups *with* inversion symmetry, you need to specify the
distribution of the electrons in the two fermion irreps
\ :sup:``[1] <#_note-saue-jensen:qsym2>`_`\ .

*Example:* Minimum wave function input for the Neon atom:

::

    **WAVE FUNCTION
    .SCF
    *SCF
    .CLOSED SHELL
     4 6

.OPEN SHELL
-----------

Specification of open shell(s).

For each open shell give the number of electrons and the number of
active spinors.

*Example:*

::

    .OPEN SHELL
     1
     5/0,6

1 open shell with 5 electrons in 6 spinors (= 3 Kramers pairs) in irrep
2. Thus, the fractional occupation is 5/6.

The open shell module in DIRAC is based on average-of-configurations
\ :sup:``[2] <#_note-thyssen-98>`_`\ . The simplest case is one electron
in two spinors (= one Kramers pair). For this special case the
average-of-configuration calculation gives the same result as the usual
restricted open-shell Hartree-Fock. For all other cases the calculation
gives the average energy of many states.

An example of one electron in two spinors is the Boron atom:

::

    **WAVE FUNCTION
    .SCF
    *SCF
    .CLOSED SHELL
     4 0
    .OPEN SHELL
     1
     1/0,2

The input is interpreted in the following way: We have 4 inactive
electrons (in the 1*s* and 2*s* orbitals) and one open shell, where 1
electron is in two spinors (the 2*p*\ :sub:`1 / 2`\  and
|2\\bar{p}\_{1/2}| spinors) in irrep 2. This calculation would give you
the *P*\ :sub:`1 / 2`\  energy of Boron. We could also do a calculation
with

::

    **WAVE FUNCTION
    .SCF
    *SCF
    .CLOSED SHELL
     4 0
    .OPEN SHELL
     1 
     1/0,6

where we have the active electrons distributed in all six *p*-orbitals.
This would give us the average energy of six states:

.. figure:: /dirac/images/math/e/7/8/e78f850c2be2f3bdd9d29321da2a7556.png
   :align: center
   :alt:  \\frac{2\\cdot E\_{P\_{1/2}} + 4\\cdot E\_{P\_{3/2}}}{6}

    \\frac{2\\cdot E\_{P\_{1/2}} + 4\\cdot E\_{P\_{3/2}}}{6}
It is also possible to specify more than one open shell. Let us take a
look at two different inputs for the open shells of the Platinum atom
(6*s*\ :sup:`1`\ 5*d*\ :sup:`9`\ ):

In the first example:

::

    .OPEN SHELL
     2
     9/10,0
     1/2,0

we get the average energy of 20 states:

|\\left(\\begin{array}{c} 10 \\\\ 9 \\end{array}\\right)\\cdot
\\left(\\begin{array}{c} 2 \\\\ 1 \\end{array}\\right)\\cdot=10\\cdot
2=20|.

All 20 states are *s*\ :sup:`1`\ *d*\ :sup:`9`\  states.

In the second input:

::

    .OPEN SHELL
     1
     10/12,0

we distribute the same ten electrons in all 12 spinors. This gives a
total of |\\left(\\begin{array}{c} 12 \\\\ 10 \\end{array}\\right) = 66|
states, which are a mixture of one *s*\ :sup:`0`\ *d*\ :sup:`10`\ 
state, 20 *s*\ :sup:`1`\ *d*\ :sup:`9`\  and 45
*s*\ :sup:`2`\ *d*\ :sup:`8`\  states.

Note that the order of closed and open shells are assumed to be as in
the following scheme:

**Virtual orbitals**

Not occupied

...

**Open shell 2**

Fractionally occupied

**Open shell 1**

Fractionally occupied

**Closed shell**

Doubly occupied

that is, the lowest molecular orbitals are doubly occupied, the next
ones are occupied with the electrons of open shell no. 1, etc.

Other orderings can be achieved by using `.REORDER
MO'S </dirac/index.php/Manual:WAVE_FUNCTION#.REORDER_MO.27S>`_ and
`.OVLSEL <#.OVLSEL>`_.

To get the energies of the individual states present in the
average-of-configurations, specify
`.RESOLVE </dirac/index.php/Manual:WAVE_FUNCTION#.RESOLVE>`_ (see also
`\*RESOLVE </dirac/index.php/Manual:WAVE_FUNCTION:RESOLVE>`_).

To get the energies of (some) of the individual states present in the
average of configurations, you can use the `GOSCIP
module </dirac/index.php/Manual:GOSCIP>`_, the `DIRRCI
module </dirac/index.php/Manual:DIRRCI>`_, or the `LUCITA
module </dirac/index.php/Manual:LUCITA>`_.

.BOSONS
-------

Occupation of boson irreps in spin-free calculation. For example, for
the D2h symmetry eight numbers in subsequent line, for the C2v symmetry
there four occupation numbers in line.

.MJSELE
-------

Occupation of MJ-splitted fermion irreps (only for the linear symmetry).
Next line contains sequence of occupation numbers, see the SCF output
for how many.

.AOC
----

Average-of-configuration calculation (default for open-shell
Hartree-Fock).

.FOCC
-----

Fractional occupation (default for open-shell
Kohn-Sham).\ :sup:`[*please `clarify </dirac/index.php/Clarify>`_*]`\ 

`.FOCC <#.FOCC>`_ calculations are less memory-intensive than
`.AOC <#.AOC>`_ calculations. In the latter case one additional AO-Fock
matrix is generated for each open shell.

`.FOCC <#.FOCC>`_ calculations are therefore an interesting option for
generating start orbitals for MCSCF as well as initial convergence in
open-shell Hartree-Fock.

.AUTOCC
-------

Program is allowed to change occupation during SCF cycles. This is
deactivated by default. However, the program will still try to do an
automatic initial occupation if neither `.CLOSED
SHELL <#.CLOSED_SHELL>`_ nor `.OPEN SHELL <#.OPEN_SHELL>`_ is given.

Trial function
==============

An SCF-calculation (HF or DFT) may be initiated in three different ways:

-  Using **MO coefficients** from a previous calculation.
-  Using coefficients obtained by diagonalization of the one-electron
   Fock matrix: the **bare nucleus approach**.
-  Using the **two-electron Fock matrix** from a previous calculation;
   this may be thought of as starting from a converged SCF potential.

The default is to start from MO coefficients if the file DFCOEF is
present. Otherwise the corrected bare nucleus approach is followed. In
all three cases linear dependencies are removed in the zeroth iteration.

.NOBNCR
-------

Switch off the bare nucleus correction. This correction is on per
default and improves the screening of the nuclei by estimating
two-electron repulsion via nuclear-attraction type integrals:

.. figure:: /dirac/images/math/6/1/b/61bb629da10a6bdd5f78bac3ad516a20.png
   :align: center
   :alt:  \\langle X\_{A} \\vert \\sum\_{C}\\frac{-Z\_{C}\\cdot \\sum\_j
   a\_j e^{(-\\alpha^C\_j r\_{C}^{2})}}{r\_{C}} \\vert X\_{B}
   \\rangle~~,~~X = L,S

    \\langle X\_{A} \\vert \\sum\_{C}\\frac{-Z\_{C}\\cdot \\sum\_j a\_j
   e^{(-\\alpha^C\_j r\_{C}^{2})}}{r\_{C}} \\vert X\_{B} \\rangle~~,~~X
   = L,S
The coefficients *a*\ :sub:`*j*`\  and the exponents |\\alpha^C\_j| in
this expression are chosen according to Slater's rules to obtain an
approximate atomic electronic density for the initial guess. For
example, with one heavy element and without this correction (that is,
with the bare nucleus Hamiltonian) all electrons will end up on that
heavy element in the initial guess!

.TRIVEC
-------

Start SCF-iterations from the vector file.

.TRIFCK
-------

Start SCF-iterations from the two-electron Fock matrix from previous
calculation (stored on file DFFCK2).

Convergence criteria
====================

Three different criteria for convergence may be chosen:

-  The norm of the DIIS error vector |\\mathbf{e} = [\\mathbf{F},
   \\mathbf{D}]| (in MO basis). This corresponds to the norm of the
   electronic gradient and is the recommended convergence criterion.
   When you are only interested in the energy `.EVCCNV <#.EVCCNV>`_ =
   1.0D-5 is usually sufficient. For properties and correlated methods
   you should converge to `.EVCCNV <#.EVCCNV>`_ = 1.0D-9. Large negative
   energy eigenvalues lead to a loss of precision that might lead to
   convergence problems. Remember also that a too loose screening
   threshold (too many integrals neglected) will hinder convergence. You
   should modify
   `.SCREEN </dirac/index.php/Manual:INTEGRALS:TWOINT#.SCREEN>`_ under
   `\*TWOINT </dirac/index.php/Manual:INTEGRALS:TWOINT>`_ if you modify
   `.EVCCNV <#.EVCCNV>`_ or one of the other two convergence criteria.
-  The difference in total energy between two consecutive iterations.
-  The largest absolute difference in the total Fock matrix between two
   consecutive iterations.

The change in total energy is approximately the square of the largest
element in the error vector or the largest change in the Fock matrix.
The default is convergence on electronic gradient with 1.0D-6 as
threshold. Alternatively, the iterations will stop at the maximum number
of iterations.

Sometimes it may happen that the specified convergence criterion is too
tight for the given basis set and/or other input parameters. In this
case one needs to decide whether one should proceed with post-HF steps
(like correlation calculations) or not. The program decides this by
looking at a secondary convergence criterion that gives the *allowed*
convergence. This value is by default the same as first or *desired*
convergence criterion but can be made lower to make sure that a
calculation does not abort when the convergence is slightly above
threshold.

For more detailed help see `SCF convergence
troubleshooting </dirac/index.php/SCF_convergence_troubleshooting>`_ or
`here </dirac/index.php/Drawbacks_and_difficult_cases_in_DIRAC#Hf2_.2B_.2C_Ta3_.2B_.2C_W4_.2B_with_broken_orbital_degeneracies>`_.

.MAXITR
-------

Maximum number of SCF iterations.

*Default:*

::

    .MAXITR
     50

When restarting SCF itrations from previous molecular orbitals file
(DFCMO or formatted DFPCMO), we recommend to decrease maximum number of
iterations together with readjusting desired and allowed convergence
thresholds. By properly set desired and allowed thresholds one can have
exact number of iterations specified by *.MAXITR*.

.EVCCNV
-------

Converge on error vector (electronic gradient).

*2 (real) Arguments:*

::

    .EVCCNV
     desired threshold allowed threshold

.ERGCNV
-------

Threshold for convergence on total energy.

*2 (real) Arguments:*

::

    .ERGCNV
     desired threshold allowed threshold

.FCKCNV
-------

Converge on largest absolute change in Fock matrix.

*2 (real) Arguments:*

::

    .FCKCNV
     desired threshold allowed threshold

Note that the *allowed* threshold may be omitted. It is then made equal
to the *desired* threshold.

Convergence acceleration
========================

It is imperative to keep the number of SCF iterations at a minimum. This
may be achieved by convergence acceleration schemes:

-  **Damping:** The simplest scheme is damping of the Fock matrix that
   may remove oscillations. In iteration *n* + 1 the Fock matrix to be
   diagonalized is: |\\mathbf{F}' = (1-c) \\mathbf{F}\_{n+1} + c
   \\mathbf{F}\_n|, where *c* is the damping factor.
-  **DIIS:** Direct inversion of iterative subspaces
   \ :sup:``[3] <#_note-pulay:diis1>`_`\ 
   \ :sup:``[4] <#_note-pulay:diis2>`_`\ 
   \ :sup:``[5] <#_note-pulay:diis3>`_`\  may be thought of as
   generalized damping involving Fock matrices from many iterations.
   Damping factors are obtained by solving a simple matrix equation
   involving the B-matrix constructed from error vectors (approximate
   gradients). Linear dependent columns in the B-matrix is removed.

In DIRAC DIIS takes precedence over damping.

.DIISTH
-------

Change the default convergence threshold for initiation of DIIS, based
on largest element of error vector.

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

.DIISMO
-------

Activate DIIS in orthogonal basis (MO) with the error vector as
described above (|image6|).

.DIISAO
-------

Activate
`DALTON <http://www.kjemi.uio.no/software/dalton/dalton.html>`_-like
DIIS using AO-basis. The error vector is |\\mathbf{e} = \\mathbf{f} -
\\mathbf{f}^{\\dagger}|

where |\\mathbf{f}| is given by

.. figure:: /dirac/images/math/4/8/d/48d4ca3b3a40eaa082bddb4ccf98e1a9.png
   :align: center
   :alt: \\mathbf{f} = \\mathbf{C}^\\dagger \\cdot
   \\mathbf{S}\_\\text{AO} \\cdot \\biggl[ \\mathbf{D}^{C}\_\\text{AO}
   \\cdot \\mathbf{F}^{D}\_\\text{AO} + \\sum\_{O \\in \\mathcal{O}}
   f\_{O} \\cdot \\mathbf{D}^{O}\_\\text{AO} \\cdot \\bigl(
   \\mathbf{F}^{D}\_\\text{AO} + ( a\_O - 1 )
   \\mathbf{Q}^{V,O}\_\\text{AO} \\bigr) \\biggr] \\mathbf{C}

   \\mathbf{f} = \\mathbf{C}^\\dagger \\cdot \\mathbf{S}\_\\text{AO}
   \\cdot \\biggl[ \\mathbf{D}^{C}\_\\text{AO} \\cdot
   \\mathbf{F}^{D}\_\\text{AO} + \\sum\_{O \\in \\mathcal{O}} f\_{O}
   \\cdot \\mathbf{D}^{O}\_\\text{AO} \\cdot \\bigl(
   \\mathbf{F}^{D}\_\\text{AO} + ( a\_O - 1 )
   \\mathbf{Q}^{V,O}\_\\text{AO} \\bigr) \\biggr] \\mathbf{C}
.NODIIS
-------

Do not perform DIIS. The default is to activate `.DIISMO <#.DIISMO>`_
for closed-shell calculations, and to activate `.DIISAO <#.DIISAO>`_ for
average-of-configurations calculations.

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
default, but DIIS takes precedence. In case all columns in the B-matrix
is removed by linear dependency, damping is activated.

Level shift
===========

.LSHIFT
-------

Activate level shift (for virtual orbitals). Followed by a real argument
(level shift).

.OLEVEL
-------

Activate level shift (for open-shell orbitals). Followed by a real
argument (level shift), one line for each open shell {Please give
example}.

State selection
===============

Convergence can be improved by selection of vectors based on overlap
with vectors from a previous iteration. This method may also be used for
convergence to some excited state.

If dynamic overlap selection is used, the vector set from the previous
iteration is used as the criterion. For the first iteration either
restart vectors or vectors generated by the bare nucleus approach
(not*recommended) are used.*

If `.NODYNSEL <#.NODYNSEL>`_ is given, either the restart vectors or the
bare nucleus vectors are used, i.e. the overlap selection vectors are
*not* updated in each iteration. Please note, that overlap selection
based on vectors from the bare nucleus approach is not recommended.

Overlap selection is very useful together with `.REORDER
MO'S </dirac/index.php/Manual:WAVE_FUNCTION#.REORDER_MO.27S>`_. This
will reshuffle the vectors within the restart coefficients. Example:
First one might do a open shell calculation on Boron, this would give
the *P*\ :sub:`1 / 2`\  state. But if we restart on the *P*\ :sub:`1 /
2`\  coefficients, interchange the *p*\ :sub:`1 / 2`\  with the
*p*\ :sub:`3 / 2`\  orbitals, and request overlap selection, we can
converge to the *P*\ :sub:`3 / 2`\  state.

There also exists a keyword for reordering the converged SCF orbitals.
This is useful for reordering the orbitals for the 4-index
transformation and subsequent correlation calculations (CCSD, CI etc.)
(see `.POST DHF REORDER
MO'S </dirac/index.php/Manual:WAVE_FUNCTION#.POST_DHF_REORDER_MO.27S>`_).

.OVLSEL
-------

Activate dynamic overlap selection. The default is no overlap selection.

.NODYNSEL
---------

No dynamic update of overlap selection vectors. The default is dynamic
update.

Iteration speedup
=================

The total run time may be reduced significantly by reducing the number
of integrals to be processed in each iteration:

-  **Screening on integrals:** Thresholds may be set to eliminate
   integrals below the threshold value
   \ :sup:``[6] <#_note-saue:csau>`_`\ . The threshold for LL integrals
   is set in the basis file, but this threshold may be adjusted for SL
   and SS integrals by threshold factors set in the
   `\*\*INTEGRALS </dirac/index.php/Manual:INTEGRALS>`_ section.

-  **Screening on density:** In direct mode further reductions are
   obtained by screening on the density matrix as well
   \ :sup:``[6] <#_note-saue:csau>`_`\ . This becomes even more
   effective if one employs *differential densities*, that is |\\Delta
   \\mathbf{D} = \\mathbf{D}\_{n+1} - \\alpha \\cdot \\mathbf{D}\_n|.
   The default value for α is |\\alpha=\\frac{ \\mathbf{D}\_{n+1} \\cdot
   \\mathbf{D}\_n }{ \\mathbf{D}\_{n} \\cdot \\mathbf{D}\_n }| which
   corresponds to a Gram-Schmidt orthogonalization. As SCF converges, α
   goes towards 1, but α can also explicitly be set equal to 1 with
   `.FIXDIF <#.FIXDIF>`_.

-  **Neglect of integrals:** The number of integrals to be processed may
   be reduced even further by adding SL and SS integrals only at an
   advanced stage in the SCF-iterations, as determined either by the
   number of iterations or by energy convergence. The latter takes
   precedence over the former.

.NODSCF
-------

Do not perform SCF-iterations with differential density matrix.

*Default:* Use differential density matrix in direct SCF.

.FIXDIF
-------

Set α equal to 1.

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
`.INTFLG </dirac/index.php/Manual:HAMILTONIAN#.INTFLG>`_ under
`\*\*HAMILTONIAN </dirac/index.php/Manual:HAMILTONIAN>`_).

*Default:* `.INTFLG </dirac/index.php/Manual:HAMILTONIAN#.INTFLG>`_ from
`\*\*HAMILTONIAN </dirac/index.php/Manual:HAMILTONIAN>`_.

Output control
==============

.PRINT
------

General print level for the SCF method. For instance, value of 2 prints
eignevalues during each iteration.

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

In studies of electronic structure it may be of interest to eliminate or
freeze certain orbitals. This option is furthermore useful for
convergence, in particular to excited electronic states. A simple case
is the thallium atom. The ground state \ :sup:`2`\ *P*\ :sub:`1 / 2`\ 
has the electronic configuration [Xe]|4f^{14}5d^{10}6s^{2}6p^1\_{1/2}|.
The first excited state \ :sup:`2`\ *P*\ :sub:`3 / 2`\  with the
electronic configuration |[Xe]4f^{14}5d^{10}6s^{2}6p^1\_{3/2}| can
easily be accessed by first calculating the ground state, then
eliminating the |6p^1\_{1/2}| from the ensuing calculation. In the final
calculation the excited state is relaxed using overlap selection. The
use of frozen orbitals is demonstrated in test/33.frozen: When the
geometry of the water molecule is optimized with the oxygen 1s and 2s
orbitals frozen, a bond angle of 96.242 degrees is found, contrary to
the 90 degrees one might have expected when s-p hybridization is thus
blocked \ :sup:``[7] <#_note-kutzelnigg:maingroup>`_`\ .

The elimination of orbitals is achieved by projecting the selected
orbitals out of the transformation matrix to orthonormal basis. The
selected orbitals can be expressed either in the full molecular basis or
in the basis set of the chosen fragment. In the latter case, the same
set of fragment orbitals can in the case of atomic fragments be used at
different molecular geometries. One may even perform a geometry
optimization, but only using the numerical gradient. When freezing
orbitals the selected orbitals are first eliminated from the
transformation to orthonormal basis, but then reintroduced in the
backtransformation step. They will appear in the output with zero
orbital eigenvalues. Note that when freezing orbitals the orbitals to be
eliminated must be specified as well. The frozen orbitals must be a
subset of the eliminated orbitals.

Fragments are defined with respect to the list of symmetry-independent
atoms appearing in the DIRAC basis file: Consider the water molecule in
the full *C*\ :sub:`2*v*`\  symmetry. Then there are two fragments: the
oxygen atom and the H\ :sub:`2`\  moiety. However, with no symmetry
there will be three fragments: the oxygen atom and the two hydrogen
atoms. At the moment there are no orthonormalization of fragments on
different fragments and so in practice one should only use orbitals from
one fragment.

.PROJECT
--------

Eliminate orbitals by projecting them out of the transformation matrix
to orthonormal basis.

*Arguments:* Number of fragments (NPRJREF).

Then for each fragment read name PRJFIL of the coefficient file followed
by the number of symmetry-independent nuclei in this fragment followed
by an `orbital string </dirac/index.php/Manual:Orbital_strings>`_ of
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

Eliminated/frozen orbitals are given in the fragment basis. Note that
the list of fragments is assumed to follow the list of symmetry
independent nuclei in the DIRAC basis file.

.FROZEN
-------

Freeze orbitals (must only be used in conjunction with
`.PROJECT <#.PROJECT>`_).

For each fermion irrep, give an `orbital
string </dirac/index.php/Manual:Orbital_strings>`_ of orbitals to
freeze.

Sample input
============

::

    *SCF
    .PRINT
     1
    .NELECT
     10 0
    .MAXITR
     30
    .EVCCNV
     1.0D-6
    .MXDIIS
     10
    .CNVINT
     1.0D-2 1.0D-4
    .EIGPRI
     1 0

References
==========

#. `↑ <#_ref-saue-jensen:qsym2_0>`_ T. Saue and H. J. Aa. Jensen,
   *Quaternion symmetry of the Dirac equation*, in *Mathematical Methods
   for Ab Initio Quantum Chemistry*, edited by M. Defrancheschi and C.
   Le Bris, Springer, Berlin, 2000.
#. `↑ <#_ref-thyssen-98_0>`_ J. Thyssen and H. J. Aa. Jensen,
   *Average-of-configurations SCF manuscript*, unpublished (1998).
#. `↑ <#_ref-pulay:diis1_0>`_ P. Pulay, *Convergence acceleration of
   iterative sequences. The case of SCF iteration*, Chem. Phys. Lett.
   **73**, *393* (1980) `electronic
   version <http://dx.doi.org/10.1016/0009-2614(80)80396-4>`_.
#. `↑ <#_ref-pulay:diis2_0>`_ P. Pulay, *Improved SCF acceleration*, J.
   Comput. Chem. **3**, *556* (1982) `electronic
   version <http://dx.doi.org/10.1002/jcc.540030413>`_.
#. `↑ <#_ref-pulay:diis3_0>`_ T. Hamilton and P. Pulay, *Direct
   inversion in the iterative subspace (DIIS) optmization of open-shell,
   excited-state, and small multiconfiguration SCF wave functions*, J.
   Chem. Phys. **84**, *5728* (1986) `electronic
   version <http://dx.doi.org/10.1063/1.449880>`_.
#. ↑ \ :sup:``6.0 <#_ref-saue:csau_0>`_`\ 
   \ :sup:``6.1 <#_ref-saue:csau_1>`_`\  T. Saue, K. Fægri, T. Helgaker,
   and O. Gropen, *Principles of direct 4-component relativistic SCF:
   Application to cesium auride*, Mol. Phys. **91**, *937* (1997)
   `electronic version <http://dx.doi.org/10.1080/002689797171058>`_.
#. `↑ <#_ref-kutzelnigg:maingroup_0>`_ W. Kutzelnigg, *Chemical Bonding
   in the Higher Main Group Elements*, Angew. Chem. Int. Ed. Engl.
   **23**, *272* (1984) `electronic
   version <http://dx.doi.org/10.1002/anie.198402721>`_.

