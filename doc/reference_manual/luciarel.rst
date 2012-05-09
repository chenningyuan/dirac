Manual:LUCIAREL
===============

**This page describes unreleased functionality. The keywords may not be
available in your version of DIRAC.**

**stefan: fixme**

*2- and 4-component relativistic GASCI module* written by Timo Fleig
based on LUCIA by Jeppe Olsen, parallelization by Stefan Knecht

`Template:LUCIAREL </dirac/index.php?title=Template:LUCIAREL&action=edit>`_
is a string-based Hamiltonian-direct configuration interaction (CI)
program, based on the LUCIA code 
[1] <#_note-olsen:lucia>
It is capable of doing efficient CI computations at arbitrary excitation
level, e.g. FCI, SDCI, RASCI, and MRCI using general active spaces. The
code is interfaced to molecular integrals obtained in the relativistic
2c- and 4-component framework and uses double-point group symmetry. It
is implemented as a full parallel version
[2] <#_note-knecht-jensen-fleig:lupar2>

A central feature of the program is the Generalized Active Space (GAS)
concept, in which the underlying total orbital space is subdivided into
a basically arbitrary number (save for an upper limit) of subspaces with
arbitrary occupation constraints. This is the most general approach to
orbital space subdivisions. The program uses DIRAC Kramers pairs from
either a closed- or an open-shell calculation in a relativistic two- and
four-componenet formalism (see
`\*\*HAMILTONIAN </dirac/index.php/Manual:HAMILTONIAN>`_).

The technical limitations are roughly set by several 100 million
determinants in the CI expansion on PCs and common computing clusters
and several billions of determinants on supercomputers with ample
memory. For calculations involving more than 500.000 determinants and in
particular higher than double excitations it is strongly recommended to
use
`Template:Luciarel </dirac/index.php?title=Template:Luciarel&action=edit>`_
in 2c- and 4c applications.

If desired, the program also computes 1- and 2-particle densities from
optimized CI wave functions. The density matrices may be printed along
with natural orbital occupations and the corresponding eigenvectors
(NOs).

Mandatory keywords
==================

.INIWFC
-------

Initial HF wave function, closed-shell:

::

    .INIWFC
     DHFSCF

or open-shell:

::

    .INIWFC
     OSHSCF

*Default:* None.

.CITYPE
-------

Type of CI calculation. Typical multi-reference CI calculations should
be defined by RASCI or GASCI.

*Values:* FCI, SDCI, SDTQ, RASCI, GASCI

*Default:* None.

.MULTIP
-------

State spin multiplicity.

*Default:* None.

Optional keywords
=================

.TITLE
------

One line with title of CI calculation.

.NROOTS
-------

Number of states to optimize on.

*Default:*

::

    .NROOTS
     1

.SYMMET
-------

State symmetry. This is the boson symmetry label referring to an irrep
ordering as defined by the group generators.

*Default:*

::

    .SYMMET
     1

.SZCALC
-------

Approximate size of calculation.

Normal CI (default):

::

    .SZCALC
     NOR

Large CI:

::

    .SZCALC
     LAR

Huge CI:

::

    .SZCALC
     HUG

.MAXITR
-------

Number of CI iterations.

*Default:*

::

    .MAXITR
     100

.INACTI
-------

Inactive orbitals per boson symmetry, separated by commas.

Not allowed with type GASCI.

*Default:* All orbitals active.

.DENSI
------

Level of computed density matrices.

For one-particle density only:

::

    .DENSI
     1

For one- and two-particle density matrices:

::

    .DENSI
     2

*Default:* calculate and print natural orbital occupation numbers.

.RSTRCI
-------

Optional restart from CI vector on file LUCVECT resp. LUCVECT.0
(**parallel calculation**):

::

    .RSTRCI
     1

*Default:* No restart.

.LBLKSZ
-------

Largest allowed vector batchsize; will be determined dynamically
otherwise.

*Default:*

::

    .LBLKSZ
     100 000

.TRUNCF
-------

Truncate residual vectors before creating new trial vector (**use with
care! experts only**).

*Default:*

::

    .TRUNCF
     1.0D-10

.PRINTL
-------

Local print level for DIRLUC interface.

*Default:*

::

    .PRINTL
     0

*Range:* 0-10

.PRINTG
-------

Global print level for LUCIA. Use with care! Set to 1 for printing
density matrices (if specified) and natural orbitals.

*Default:*

::

    .PRINTG
     0

*Range:* 0-4

Specific input -- GASCI
=======================

.NACTEL
-------

Number of active electrons.

*Default:* None.

.GASSHE
-------

Number and specification of GAS orbitals. Line with the number of GA
spaces used (1-16), followed by one line per GAS with number of orbitals
per boson symmetry, separated by commas.

*Default:* None.

.GASSPC
-------

Number and specification of sequential CI calculations. Line with the
number of CI calculations with given GA spaces (currently only 1 is
allowed), followed by one line per GAS with 2 numbers each: The first
gives the minimal number of accumulated electrons after this GAS, the
second the corresponding maximum number, separated by blanks (defining
occupation constraints of each GAS).

*Default:* None.

Specific input -- RASCI
=======================

.NACTEL
-------

Number of active electrons.

*Default:* None.

.INACTI
-------

Inactive orbitals per boson symmetry, separated by commas.

*Default:* All orbitals active.

.RAS1
-----

RAS1 specification and maximum number of holes. Line with orbitals per
boson symmetry, separated by commas, followed by a line with the maximum
number of holes in RAS1.

*Default:* None.

.RAS2
-----

RAS2 specification. Line with orbitals per boson symmetry, separated by
commas.

*Default:* None.

.RAS3
-----

RAS3 specification and maximum number of electrons. Line with orbitals
per boson symmetry, separated by commas, followed by a line with the
maximum number of electrons in RAS3.

*Default:* None.

Parallel keywords
=================

.DISTRT
-------

Specifies vector block distribution routine to use in parallel
calculation.

simple distribution routine useful in small calculations:

::

    .DISTRT
     1

more advanced and efficient distribution routine:

::

    .DISTRT
     2

*Default:* 2

References
==========

<#_ref-olsen:lucia_0>`_ J. Olsen and P. Joergensen and J. Simons,
*Passing the one-billion limit in Full CI calculations*, Chem. Phys.
Lett. **169**, *463* (1990) .

<#_ref-knecht-jensen-fleig:lupar2_0>`_ S. Knecht and H. J. Aa.
Jensen and T. Fleig, *Large-Scale Parallel Configuration Interaction.
II. Two- and 4-Component Double-Group General Active Space
Implementation with Application to BiH.*, J. Chem. Phys. **128**,
*014108* (2010) `electronic
version <http://dx.doi.org/http://link.aip.org/link/JCPSA6/v132/i1/p014108/s1>`_.

