

.. _integrals:
.. index:: **INTEGRALS

=============
\*\*INTEGRALS
=============


General directives
==================


.. index:: .NUCMOD
.NUCMOD
-------

Specify nuclear model.

Point nucleus::

  .NUCMOD
   1

Gaussian charge distribution (default)::

  .NUCMOD
   2

The point nucleus model is useful to compare Lévy-Leblond type calculations
with regular nonrelativistic calculations done with another code, e.g. `Dalton
<http://www.daltonprogram.org>`_.  The two methods should give precisely the
same energies.

For the Gaussian charge distribution the default exponents are in
accordance with values proposed by Visscher and Dyall :cite:`Visscher1997b`.


.. index:: .SELECT
.SELECT
-------

Restrict range of nuclei in one-electron integrals involving single
atomic centers, for example electric field gradients.

Example: Restrict range to four nuclei (first number), the nuclei 1, 3, 7, and 8::

  .SELECT
   4
   1
   3
   7
   8


.. index:: .PHASEORIGIN
.PHASEORIGIN
------------

Origin appearing in the London atomic orbital phase-factors. Default::

  .PHASEORIGIN
   0.0 0.0 0.0

.. warning::

  This is the text from the DALTON manual. We wonder if it is used for
  other integrals, in the code in her1int.F it is the final "else" option.


.. index:: .GAUGEORIGIN
.GAUGEORIGIN
------------

Set the gauge origin. By default it is set to the coordinate origin::

  .GAUGEORIGIN
   0.0 0.0 0.0


.. index:: .DIPORG
.DIPORG
-------

Origin for all moment integrals, including the dipole (dipole is independent of
origin only for neutral systems). By default it is set to the coordinate
origin::

  .DIPORG
   0.0 0.0 0.0


.. index:: .MAGCOR
.MAGCOR
-------

Print the symmetrized nuclear magnetic moments. This corresponds to taking
symmetry combinations of rotations, not coordinates, at each nuclear center.
The numbering is used in labels of various magnetic integrals.


.. index:: .PRINT
.PRINT
------

General print level in the integral module. Default::

  .PRINT
   1


.. index:: *ONEINT
\*ONEINT
========

This subsection gives directives for the generation of one-electron
integrals. Based on input in the other section, the program will
determine what integrals to calculate.


.. index:: .PRINT
.PRINT
------

Print level in one-electron integral routines. By default
print level is taken from \*\*INTEGRALS.


.. index:: *READIN
\*READIN
========

This subsection allows changes of defaults in the reading of the mol file.


.. index:: .UNCONTRACT
.UNCONTRACT
-----------

Decontract basis sets specified as contracted (when being read from the
library, for example). In the case of using non-relativistically contracted
basis set the decontraction is necessary for heavy elements.  Two-component
quasirelativistic Hamiltonians (like X2C) work only with decontracted basis
set. By default only the small component is decontracted.


.. index:: .PRINT
.PRINT
------

Print level in the reading of the mol file. By default
print level is taken from \*\*INTEGRALS.


.. index:: .MAXPRI
.MAXPRI
-------

Maximum number of primitive functions in a given block in basis file.
Default::

  .MAXPRI
   22


.. index:: *TWOINT
\*TWOINT
========

This subsection gives directives for the generation of two-electron integrals.
It also gives directives for the construction of Fock matrices, such as
screening.


.. index:: .SCREEN
.SCREEN
-------

Screening threshold for integral direct calculations of Fock matrices :cite:`Saue1997`.
Default::

  .SCREEN
   1.0D-12

Note that the screening threshold may influence the convergence. In general,
the screening threshold must be about three orders of magnitude smaller than
the desired norm of the electronic gradient at convergence.

Choosing a negative value for the screening turns the integral screening
completely off.


.. index:: .ICEDIF
.ICEDIF
-------

Separate screening of Coulomb and exchange contributions :cite:`Saue1997`.
Useful for fine regulation of the convergence process.
Default: Coulomb and exchange on (1 = on; off = 0)::

  .ICEDIF
   1 1


.. index:: .THRFAC
.THRFAC
-------

Adjust the integral thresholds for SL and SS integrals. For conventional
integral calculations only integrals above the threshold given in the mol file
are written to disk. The thresholds for the SL and SS integrals are divided by
the factors given here. Default::

  .THRFAC
   1.0 1.0


.. index:: .SOFOCK
.SOFOCK
-------

Do direct Fock matrix construction in symmetry-adapted basis (SO basis).

By default the direct Fock matrix construction in performed AO basis using the
skeleton matrix approach. This may give better screening, but is more memory
intensive.


.. index:: .PRINT
.PRINT
------

Set the print level in two-electron integral routines for the
calculation of a particular shell quadruplet. The print level is changed
only for the given shell quadruplet. A zero matches all shells, thus::

  .PRINT
   4 0 0 0 0

or just::

  .PRINT
   4

sets the print level to 4 for all shell quadruplets.

Use with care to avoid massive output! At print level 15 the individual
integrals are printed.


.. index:: .TIME
.TIME
-----

Give detailed timing for integral calculation.
