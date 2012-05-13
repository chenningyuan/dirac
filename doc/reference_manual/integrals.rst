\*\*INTEGRALS -- Integral directives
====================================

.NUCMOD
-------

Nuclear model.

*Point nucleus:*

::

    .NUCMOD
     1

*Gaussian charge distribution (default):*

::

    .NUCMOD
     2

The point nucleus model is useful to compare Lévy-Leblond type
calculations with regular nonrelativistic calculations done with another
code, e.g.
`DALTON <http://www.kjemi.uio.no/software/dalton/dalton.html>`_. The two
methods should give precisely the same energies.

For the Gaussian charge distribution the default exponents are in
accordance with values proposed by Visscher and Dyall [Visscher1997b]_
see also the web-site
`[1] <http://dirac.chem.sdu.dk/doc/FiniteNuclei/FiniteNuclei.shtml>`_.
If desired, other exponents may be specified in the `basis
file </dirac/index.php/Manual:Basis:Title_and_symmetry>`_.

**Advanced options**

.SELECT
-------

Restrict range of nuclei in one-electron integrals involving single
atomic centers, for example electric field gradients.

*Example:* Restrict range to four nuclei, the nuclei 1, 3, 7, and 8.

::

    .SELECT
     4
     1
     3
     7
     8

.PHASEO
-------

Origin appearing in the London atomic orbital phase-factors.

*Default:*

::

    .PHASEO
     0.0 0.0 0.0

This is the text from the DALTON manual. We wonder if it is used for
other integrals, in the code in her1int.F it is the final "else" option.

.GAUGEO
-------

Gauge origin.

*Default:*

::

    .GAUGEO
     0.0 0.0 0.0

.DIPORG
-------

Origin for all moment integrals, including the dipole (dipole is
independent of origin only for neutral systems).

*Default:*

::

    .DIPORG
     0.0 0.0 0.0

.MAGCOR
-------

Print the symmetrized nuclear magnetic moments; this corresponds to
taking symmetry combinations of rotations, not coordinates, at each
nuclear center. The numbering is used in labels of various magnetic
integrals.

**Programmers options**

.PRINT
------

General print level in HERMIT.

*Default:*

::

    .PRINT
     1

\*ONEINT -- One-electron integrals
==================================

This subsection gives directives for the generation of one-electron
integrals. Based on input in the other section, the program will
determine what integrals to calculate.

**Programmers options**

.PRINT
------

Print level in one-electron integral routines (abacus).

*Default:* Print level taken from
`\*\*INTEGRALS </dirac/index.php/Manual:INTEGRALS>`_.


Manual:INTEGRALS:READIN
=======================

\*READIN -- Reading basis set file
==================================

This subsection allows changes of defaults in the reading of the `basis
file </dirac/index.php/Manual:Basis:Title_and_symmetry>`_.

**Advanced options**

.UNCONTRACT
-----------

Decontract basis sets specified as contracted (when being read from the
library, for example). In the case of using non-relativistically
contracted basis set the decontraction is necessary for heavy elements.
Two-component quasirelativistic Hamiltonians (like IOTC) work only with
decontracted basis set.

*Default:* Decontract only the small component.

**Programmers options**

.PRINT
------

Print level in the reading of the basis file.

*Default:* Print level from
`\*\*INTEGRALS </dirac/index.php/Manual:INTEGRALS>`_.

.MAXPRI
-------

Maximum number of primitive functions in a given block in basis file.

*Default:*

::

    .MAXPRI
     22


Manual:INTEGRALS:TWOINT
=======================

\*TWOINT -- Two-electron integrals
==================================

This subsection gives directives for the generation of two-electron
integrals.It also gives directives for the construction of Fock
matrices, such as screening.

**Advanced options**

.SCREEN
-------

Screening threshold for integral direct calculations of Fock matrices
[1] <#_note-saue:csau>

*Default:*

::

    .SCREEN
     1.0D-12

Note that the screening threshold may influence the convergence. In
general, the screening threshold must be about three orders of magnitude
smaller than the desired norm of the electronic gradient at convergence.

Choosing a negative value for the screening turns the integral screening
completely off.

.ICEDIF
-------

Separate screening of Coulomb and exchange contributions [Saue1997]_.
Useful for fine regulation of the
convergence process.

*Default:* Coulomb and exchange on (1 = on; off = 0).

::

    .ICEDIF
     1 1

.THRFAC
-------

Adjust the integral thresholds for SL and SS integrals. For conventional
integral calculations only integrals above the threshold given in the
`basis file </dirac/index.php/Manual:Basis:Title_and_symmetry>`_ are
written to disk. The thresholds for the SL and SS integrals are divided
by the factors given here.

*Default:*

::

    .THRFAC
     1.0 1.0

.SOFOCK
-------

Do direct Fock matrix construction in symmetry-adapted basis (SO basis).

*Default:* Direct Fock matrix construction in AO basis using the
skeleton matrix approach. This may give better screening, but is more
memory intensive.

**Programmers options**

.PRINT
------

Set the print level in two-electron integral routines for the
calculation of a particular shell quadruplet. The print level is changed
only for the given shell quadruplet. A zero matches all shells, thus

::

    .PRINT
     4 0 0 0 0

or just

::

    .PRINT
     4

sets the print level to 4 for all shell quadruplets.

Use with care to avoid massive output! At print level 15 the individual
integrals are printed.

.TIME
-----

Give detailed timing for integral calculation.


Manual:INTEGRALS:Sample inputs
==============================

Use the point nucleus model to compare with nonrelativistic programs
(most of them use point nucleus):

::

    **DIRAC
    .TITLE
      INTEGRALS sample input 1
    .WAVE FUNCTION
    **INTEGRALS
    .NUCMOD
     1
     [...]

Read very large contracted basis set (from a basis set file, for
example) as uncontracted.

::

    **DIRAC
    .TITLE
      INTEGRALS sample input 2
    .WAVE FUNCTION
    **INTEGRALS
    *READIN
    .UNCONTRACT
    ! comment: for demonstration we have 30s25p18d/22s18p12d...
    ! ...consumes a lot of memory
    .MAXPRI
     30 
     [...]

