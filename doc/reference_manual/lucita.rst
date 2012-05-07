Manual:LUCITA
=============

*Spin-free relativistic GASCI module* written by Jeppe Olsen, adaptation
to DIRAC by Timo Fleig, parallelization by Stefan Knecht

LUCITA is a string-based Hamiltonian-direct configuration interaction
(CI) program, based on the LUCIA code
[1] <#_note-olsen:lucia>`_`\ . 
It is capable of doing efficient
CI computations at arbitrary excitation level, e.g. FCI, SDCI, RASCI,
and MRCI using general active spaces. The code is interfaced
[2] <#_note-fleig-visscher:au2>`_`\ 
to molecular integrals
obtained in the spin-free Dirac formalism and uses non-relativistic
point group symmetry. It is implemented as a full parallel version
[3] <#_note-knecht-jensen-fleig:lupar1>`_`\ .

A central feature of the program is the Generalized Active Space (GAS)
concept, in which the underlying total orbital space is subdivided into
a basically arbitrary number (save for an upper limit) of subspaces with
arbitrary occupation constraints. This is the most general approach to
orbital space subdivisions. The program uses DIRAC orbitals from either
a closed- or an open-shell calculation in a spin-free relativistic
formalism or the non-relativistic Lévy-Leblond formalism (see
`\*\*HAMILTONIAN </dirac/index.php/Manual:HAMILTONIAN>`_).

The technical limitations are roughly set by several 100 million
determinants in the CI expansion on PCs and common computing clusters
and several billions of determinants on supercomputers with ample
memory. For calculations involving more than 500.000 determinants and in
particular higher than double excitations it is strongly recommended to
use LUCITA in spin-free DIRAC applications.

If desired, the program also computes 1- and 2-particle densities from
optimized CI wave functions. The density matrices may be printed along
with natural orbital occupations and the corresponding eigenvectors
(NOs).

LUCITA is provided with an interface section DIRLUC incorporated in the
code. The DIRLUC environment reads the DIRAC input and converts it to a
LUCIA input which is processed through an ASCII file (LUCIAIN). The
types of CI which may be specified by keyword are: Full CI (FCI), CI
with single and double excitations (SDCI), CI with single, double,
triple, and quadruple excitations (SDTQ), restricted active space CI
(RASCI), and generalized active space CI (GASCI). All types of CI can be
considered special cases of a GASCI, so internally, DIRLUC transfers all
information into GAS format. More detailed information on GAS-type of
wavefunctions can be found in the literature
<#_note-fleig-olsen-visscher:luciarel>`_`\ .


*Important Information for Users*
---------------------------------

In DIRAC08 the default integral transformation scheme has changed.
Please add the lines

::

    .SCHEME
     4

to your `MOLTRA </dirac/index.php/Manual:MOLTRA>`_ input section for the
integral transformation.

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

<#_ref-fleig-visscher:au2_0>`_ T. Fleig and L. Visscher,
*Large-Scale Electron Correlation Calculations in the Framework of
the Spin-Free Dirac Formalism. The Au\ :sub:`2`\  Molecule
Revisited*, Chem. Phys. **311**, *113* (2005) `electronic
version <http://dx.doi.org/10.1016/j.chemphys.2004.10.003>`_.

<#_ref-knecht-jensen-fleig:lupar1_0>`_ S. Knecht and H. J. Aa.
Jensen and T. Fleig, *Large-Scale Parallel Configuration Interaction.
I. Non-Relativistic and Scalar-Relativistic General Active Space
Implementation with Application to* (Rb-Ba)\ :sup:`+`\ .*, J. Chem.
Phys. **128**,*014108*(2008)* `electronic
version <http://dx.doi.org/10.1063/1.2805369>`_.

<#_ref-fleig-olsen-visscher:luciarel_0>`_ T. Fleig and J. Olsen
and L. Visscher, *The generalized active space concept for the
relativistic treatment of electron correlation: II. Large-scale
configuration interaction implementation based on relativistic 2-and
4-spinors and its application.*, J. Chem. Phys. **119**, *2963*
(2003) `electronic version <http://dx.doi.org/10.1063/1.1590636>`_.



Manual:LUCITA:Sample inputs
===========================

SDCI:

::

    *LUCITA
    .INIWFC
      DHFSCF
    .CITYPE
      SDCI
    .MULTIP
      3
    *END OF

GASCI:

::

    *LUCITA
    .TITLE
      Aluminum atom in C2h symmetry, 7s5p2d (L) basis, SDT/SD CI (1s inactive).
    .INIWFC
      OSHSCF
    .CITYPE
      GASCI
    .MULTIP
      2
    .NACTEL
      13
    .GASSHE
     5
     1,0,0,0    ! 1s
     2,0,0,0    ! 2s3s
     0,0,2,4    ! 2p3p
     5,2,1,2    ! 4s5s 4p 3d
     5,2,2,4    ! 6s7s 5p6p 4d
    .GASSPC
     1
      2  2
      3  6
     11 13
     11 13
     13 13
    *END OF

