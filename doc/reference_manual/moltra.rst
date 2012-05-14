\*\*MOLTRA -- Integral transformation module
============================================

This module transforms integrals from the scalar AO basis to integrals
in the molecular spinor basis. The algorithm is not yet described in
detail but it may be helpful to read Ref. [Visscher2002]_
to find information
about the double quaternion formalism that is used. This module is
usually invoked by one of the correlation modules but can also be run as
a stand-alone module by specifying
`.4INDEX </dirac/index.php/Manual:WAVE_FUNCTION#.4INDEX>`_ in the
`\*\*DIRAC </dirac/index.php/Manual:DIRAC>`_ section of input. The input
serves to change the default active space and to specify special wishes
(e.g modification of `.INTFLG <#.INTFLG>`_).

Three separate transformations are considered:

-  Transformation of a core Fock matrix, giving effective one-electron
   matrices.
-  Transformation of two-electron Coulomb integrals.
-  Transformation of integrals over general (one-electron) operators.

The output is written to the files MDCINT, MRCONEE, and MDPROP,
respectively, and can be read in by the
`RELCCSD </dirac/index.php/Manual:RELCCSD>`_ and
`DIRRCI </dirac/index.php/Manual:DIRRCI>`_ programs to perform
correlated calculations and by the
`RELADC </dirac/index.php/Manual:RELADC>`_ program to perform propagator
calculations. Some common input needs to be specified both in this
section and in the `DIRRCI </dirac/index.php/Manual:DIRRCI>`_ sections.

.ACTIVE
-------

Specify the active set of spinors.

*Arguments:* One line with an `orbital
string </dirac/index.php/Manual:Orbital_strings>`_ for each fermion
irrep.

*Default:*

::

    energy -20.0 10.0 1.0

**Advanced options**

.THROUT
-------

Threshold value for writing transformed 2-electron integrals to file. A
higher value will decrease the size of the MDCINT files but will
decrease the accuracy of the calculation.

*Default:*

::

    .THROUT
     1.0D-14

.PRPTRA
-------

Perform transformation of property integrals. The type(s) of integrals
that are to be transformed can be specified with the .OPERATOR keyword.

*Default:* Do not transform property integrals.

.SCHEME
-------

Selects the integral transformation algorithm ("scheme"). Several
schemes have been implemented, but users are to consider only "scheme 4"
and "scheme 6" for production calculations. Other schemes refer to
outdated to special-purpose algorithms, and should be used with care.
The most suitable scheme to use depends to a significant extent on the
type of calculation in question (size of the active space etc), and on
the computational resources used (compilers, amount of memory/scratch
disk space, whether or not disks are shared to local to compute nodes
etc).

*Scheme 6 (Default):*

::

    .SCHEME
     6

Scheme 6 is designed to decrease I/O operations at the expense of more
communication between processes, something which is more often than not
advantageous in today's architectures. It should be noted, however, that
the significant amount of communication between processes happens for
the so-called "1HT" step ("1HT" standing for first half-transformation,
since it is where the first two indexes of the electron repulsion
integrals are transformed to MO basis and these intermediate,
half-transformed entities are stored on disk), with no communication on
the "2HT" step (which will yield the final integrals in MO basis). Due
to the large amount of half-transformed entities, this scheme may demand
large scratch spaces (e.g. a Dirac-Coulomb calculation on water with the
aug-cc-pV5Z basis, where the full spinor space is active and LL and SL
integrals are transformed can reach nearly 1Tb of disk use) but with the
advantage that the full size of the problem is exactly divided between
processes - so the larger the number of MPI processes, the less disk
space is used per process.

**(Not part of the Dirac08 release)** In connection to scheme 6 the
user can also choose whether or perform an intermediate sorting step
prior to the second half-transformation step, via the keyword
**.HTSORT**. In this step the 1HT entities on disk are rearranged from
their "natural" order (arising from the two-electron integral algorithm)
to an order which will minimize the amount of I/O operations taking
place at the 2HT step.

 *HT sort disabled (Default):*

::

    .SCHEME
     6
    .HTSORT
     0

**(Not part of the Dirac08 release)** Sorting of HT intermediates is
disabled by default. Without such sorting, scheme 6's performance can be
negatively affected on computer systems where reading data is relatively
expensive (one example is for IBM p4/p5/p6 systems and the XLF compiler)
even for relatively modest calculations Dirac-Coulomb calculations
(where LL and SL integrals are considered). However, in two-component
calculations (where only LL integrals are transformed and therefore
there are much less 1HT entities on disk) of similar or larger size, in
terms of active space, good performance can still be achieved in the
same architecture.

 *HT sort (strategy #1) enabled:*

::

    .SCHEME
     6
    .HTSORT
     1

**(Not part of the Dirac08 release)** This sort strategy minimizes the
amount of read/write operations by taking up as much 1HT data in memory
and sorting it before writing it to disk again. This provides efficiency
for systems (such as the IBM systems mentioned above) where reading is
expensive, at the expense of larger memory usage (rule of thumb: the
extra memory required \*per MPI process\* will be about 1-2% of the
\*total disk\* space taken up by the 1HT quantities - so if that is 1Tb,
the extra memory used will be of of about 1-2Gb; see the calculation
output for more details).

 *Scheme 4:*

::

    .SCHEME
     4

Scheme 4 is designed to reduce communication between different
processors at the expense of increasing I/O operations. This increase in
I/O operations stems from the integrals being replicated for all MPI
processes, so that doubling the number of MPI processes doubles the
total disk usage but not the usage per MPI process.

.INTFLG
-------

Specify which classes of integrals should contribute to the transformed
two-electron integrals and the effective (core) Hamiltonian. Default is
to follow the definitions given in \*\*HAMILTONIAN. If you want
different integral classes in the 2- and 4-index transformation, then
use the `.INTFL2 <#.INTFL2>`_ and `.INTFL4 <#.INTFL4>`_ keywords below
(see `.INTFLG </dirac/index.php/Manual:HAMILTONIAN#.INTFLG>`_ under
`\*\*HAMILTONIAN </dirac/index.php/Manual:HAMILTONIAN>`_).

*Default:* `.INTFLG </dirac/index.php/Manual:HAMILTONIAN#.INTFLG>`_ from
`\*\*HAMILTONIAN </dirac/index.php/Manual:HAMILTONIAN>`_.

.INTFL2
-------

Specify what two-electron integrals to include in the 2-index
transformation ,i.e., in the two-electron part of the Fock core matrix
(see `.INTFLG </dirac/index.php/Manual:HAMILTONIAN#.INTFLG>`_ under
`\*\*HAMILTONIAN </dirac/index.php/Manual:HAMILTONIAN>`_).

*Default:* `.INTFLG </dirac/index.php/Manual:HAMILTONIAN#.INTFLG>`_ from
`\*\*HAMILTONIAN </dirac/index.php/Manual:HAMILTONIAN>`_.

.INTFL4
-------

Specify what two-electron integrals to include in the 4-index
transformation (see
`.INTFLG </dirac/index.php/Manual:HAMILTONIAN#.INTFLG>`_ under
`\*\*HAMILTONIAN </dirac/index.php/Manual:HAMILTONIAN>`_).

*Default:* `.INTFLG </dirac/index.php/Manual:HAMILTONIAN#.INTFLG>`_ from
`\*\*HAMILTONIAN </dirac/index.php/Manual:HAMILTONIAN>`_.

.CORE
-----

Specify the frozen core spinors,

For each fermion irrep, give an `orbital
string </dirac/index.php/Manual:Orbital_strings>`_ of active orbitals.

*Default:* All occupied orbitals that are not active.

.POSITRONS
----------

Allow negative energy ("positronic") spinors in active set. If this
keyword is not specified then only electronic spinors will be included
and any positronic spinors specified with `.ACTIVE <#.ACTIVE>`_ are
ignored.

.NO2IND
-------

Skip the 2-index transformation of the effective Fock matrix.

.NO4IND
-------

Skip the 4-index transformation.

.SCREEN
-------

Screening threshold in 4-index transformation (a negative value disables
screening).

*Default:*

::

    .SCREEN
     1.0D-14

.ASCII
------

(feature will become available in DIRAC10)

Write integrals to the ASCII file MO\_integrals.txt to facilitate
interfaces with correlation codes. Warning: only meant for initial
interfacing as this format is very inefficient.

*Default:* Do not write the ASCII file.

.RCORBS
-------

Recanonize orbitals before transforming.

*Default:* Do not recanonize orbitals before transforming.

**Programmers options**

.PRINT
------

Print level.

*Default:*

::

    .PRINT
     0

.2INDEX
-------

Ranges of active orbitals in 2-index transformation module specified for
index 1 and 2 and fermion irrep 1 and 2.

*Default:* Ranges set by `.ACTIVE <#.ACTIVE>`_.

.4INDEX
-------

Ranges of active orbitals in 4-index transformation module specified for
index 1 to 4 and fermion irrep 1 and 2.

*Default:* Ranges set by `.ACTIVE <#.ACTIVE>`_.

.MDCINT
-------

Write 4-index transformed integrals to file (default).

.NOMDCI
-------

Skip writing 4-index transformed integrals to file. Instead, write file
MDINFO, which contains infos used in the MCSCF program.

.SCATTER
--------

Scatter 4-index transformed integrals to all nodes (default).

.NOSCAT
-------

Do not scatter 4-index transformed integrals to all nodes

.PAR4BAS
--------

Specify number of batches for each node in 4-index transformation. Only
working for parallel distributions.

*Default:*

::

    .PAR4BAS
     -1

\*PRPTRA -- Property integrals transformation
=============================================

This section defines the property operators that should to be
transformed. Such integrals are only used in an experimental section of
`RELCCSD </dirac/index.php/Manual:RELCCSD>`_, and this subsection can
thus usually be omitted.

**Programmers options**

.PRINT
------

Print level.

*Default:*

::

    .PRINT
     0

.OPERATOR
---------

Specification of general one-electron operators.

See the `specification of one-electron
operators </dirac/index.php/Manual:One-electron_operators>`_ section for
more information and explicit examples.

