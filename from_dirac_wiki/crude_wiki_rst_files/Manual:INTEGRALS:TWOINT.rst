Manual:INTEGRALS:TWOINT
=======================

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

Contents
--------

-  `1 \*TWOINT -- Two-electron
   integrals <#.2ATWOINT_--_Two-electron_integrals>`_

   -  `1.1 .SCREEN <#.SCREEN>`_
   -  `1.2 .ICEDIF <#.ICEDIF>`_
   -  `1.3 .THRFAC <#.THRFAC>`_
   -  `1.4 .SOFOCK <#.SOFOCK>`_
   -  `1.5 .PRINT <#.PRINT>`_
   -  `1.6 .TIME <#.TIME>`_
   -  `1.7 References <#References>`_

\*TWOINT -- Two-electron integrals
==================================

This subsection gives directives for the generation of two-electron
integrals.It also gives directives for the construction of Fock
matrices, such as screening.

`|image0| </dirac/index.php/Image:Dirac.png>`_

**Advanced options**

.SCREEN
-------

Screening threshold for integral direct calculations of Fock matrices
\ :sup:``[1] <#_note-saue:csau>`_`\ .

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

Separate screening of Coulomb and exchange contributions
\ :sup:``[1] <#_note-saue:csau>`_`\ . Useful for fine regulation of the
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

`|image1| </dirac/index.php/Image:HJAaJ.png>`_

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

References
----------

#. ↑ \ :sup:``1.0 <#_ref-saue:csau_0>`_`\ 
   \ :sup:``1.1 <#_ref-saue:csau_1>`_`\  T. Saue, K. Fægri, T. Helgaker,
   and O. Gropen, *Principles of direct 4-component relativistic SCF:
   Application to cesium auride*, Mol. Phys. **91**, *937* (1997)
   `electronic version <http://dx.doi.org/10.1080/002689797171058>`_.

