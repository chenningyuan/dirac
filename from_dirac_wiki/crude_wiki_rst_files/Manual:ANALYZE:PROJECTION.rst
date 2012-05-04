Manual:ANALYZE:PROJECTION
=========================

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

Contents
--------

-  `1 \*PROJECTION -- Projection
   analysis <#.2APROJECTION_--_Projection_analysis>`_

   -  `1.1 .VECPRJ <#.VECPRJ>`_
   -  `1.2 .VECREF <#.VECREF>`_
   -  `1.3 .OWNBAS <#.OWNBAS>`_
   -  `1.4 .PROTHR <#.PROTHR>`_
   -  `1.5 .WGPOP <#.WGPOP>`_
   -  `1.6 .PRINT <#.PRINT>`_
   -  `1.7 References <#References>`_

\*PROJECTION -- Projection analysis
===================================

The current solution may be projected down onto another set of
coefficients generated from the basis, e.g. one may project molecular
solutions down onto atomic solutions in order to evaluate atomic
contributions. The coefficients of the fragments are read simultaneously
and the overlap between them is taken into account
\ :sup:``[1] <#_note-saue:TlAt>`_`\ .

The projection analysis can be thought of as a Mulliken analysis based
on the atomic (fragment) orbitals; it is therefore very much less
sensitive to the basis set.

Normally, the fragments are calculated in the full molecular basis. This
is done by zeroing charges of the remaining atoms of the molecule in the
basis set input and adjusting occupation in the menu file. It is,
however, possible to calculate the fragments in the own basis (a subset
of the full molecular basis) using `.OWNBAS <#.OWNBAS>`_. The advantage
is faster calculations and conservation of atomic symmetry for atomic
fragments. To avoid working with symmetry-combinations of atomic centers
for the fragments, it may sometimes be advantageous to dump molecular
coefficients in *C*\ :sub:`1`\  symmetry using
`.ACMOUT </dirac/index.php/Manual:GENERAL#.ACMOUT>`_ under
`\*\*GENERAL </dirac/index.php/Manual:GENERAL>`_ and do the analysis
without symmetry. When `.OWNBAS <#.OWNBAS>`_ is used, one then needs to
calculate each atomic type only once in its own basis in order to do the
complete analysis.

.VECPRJ
-------

For each fermion irrep, give an `orbital
string </dirac/index.php/Manual:Orbital_strings>`_ of orbitals to
analyze.

*Default:* Analyze the occupied electronic solutions.

.VECREF
-------

First give number of fragments to project onto, and then for each
fragment give filename of MO coefficients and the number of
symmetry-independent nuclei in this fragment and for each fermion irrep,
give an `orbital string </dirac/index.php/Manual:Orbital_strings>`_ of
reference orbitals.

*Example:*

::

    .VECREF
    4
    H00001
    1
    1
    H00002
    1
    1
    X00003
    1
    1..43
    X00004
    1
    1..43

.OWNBAS
-------

Calculate fragments in their own basis.

This keyword must be used with some care as the list of fragments is
assumed to be identical to that of symmetry independent centers.

`|image0| </dirac/index.php/Image:Dirac.png>`_

**Advanced options**

.PROTHR
-------

Set threshold for absolute value of projection coefficients to be
printed.

*Default:*

::

    .PROTHR
     0.01

.WGPOP
------

Split overlap densities according to weight of contributions.

`|image1| </dirac/index.php/Image:HJAaJ.png>`_

**Programmers options**

.PRINT
------

Print level.

*Default:*

::

    .PRINT
     0

References
----------

#. `↑ <#_ref-saue:TlAt_0>`_ K. Fægri and T. Saue, *Diatomic molecules
   between very heavy elements of group 13 and group 17: A study of
   relativistic effects on bonding*, J. Chem. Phys. **115**, *2456*
   (2001) `electronic version <http://dx.doi.org/10.1063/1.1385366>`_.

