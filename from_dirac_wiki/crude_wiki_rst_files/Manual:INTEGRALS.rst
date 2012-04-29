Manual:INTEGRALS
================

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

Contents
--------

-  `1 \*\*INTEGRALS -- Integral
   directives <#.2A.2AINTEGRALS_--_Integral_directives>`_

   -  `1.1 .NUCMOD <#.NUCMOD>`_
   -  `1.2 .SELECT <#.SELECT>`_
   -  `1.3 .PHASEO <#.PHASEO>`_
   -  `1.4 .GAUGEO <#.GAUGEO>`_
   -  `1.5 .DIPORG <#.DIPORG>`_
   -  `1.6 .MAGCOR <#.MAGCOR>`_
   -  `1.7 .PRINT <#.PRINT>`_
   -  `1.8 References <#References>`_

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
accordance with values proposed by Visscher and Dyall
\ :sup:``[1] <#_note-visscher-dyall:nucmod>`_`\ , see also the web-site
`[1] <http://dirac.chem.sdu.dk/doc/FiniteNuclei/FiniteNuclei.shtml>`_.
If desired, other exponents may be specified in the `basis
file </dirac/index.php/Manual:Basis:Title_and_symmetry>`_.

`|image0| </dirac/index.php/Image:Dirac.png>`_

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

`|image1| </dirac/index.php/Image:HJAaJ.png>`_

**Programmers options**

.PRINT
------

General print level in HERMIT.

*Default:*

::

    .PRINT
     1

References
----------

#. `↑ <#_ref-visscher-dyall:nucmod_0>`_ L. Visscher and K. G. Dyall,
   *Dirac-Fock atomic electronic structure calculations using different
   nuclear charge distributions*, At. Data Nucl. Data Tabl. **67**,
   *207* (1997) `electronic
   version <http://dx.doi.org/10.1006/adnd.1997.0751>`_.

