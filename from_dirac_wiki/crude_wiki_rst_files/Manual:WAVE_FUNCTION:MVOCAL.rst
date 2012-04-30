Manual:WAVE FUNCTION:MVOCAL
===========================

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

Contents
--------

-  `1 \*MVOCAL -- Modified virtual
   orbitals <#.2AMVOCAL_--_Modified_virtual_orbitals>`_

   -  `1.1 .PRINT <#.PRINT>`_
   -  `1.2 .VECMVO <#.VECMVO>`_
   -  `1.3 .INTFLG <#.INTFLG>`_
   -  `1.4 References <#References>`_

\*MVOCAL -- Modified virtual orbitals
=====================================

The virtual canonical Hartree-Fock orbitals are not the optimal choice
for correlated calculations in which the virtual space is truncated; the
lower orbitals tend to be diffuse since they see the *N* + 1 electron
system. One may exploit the freedom of rotation amongst the virtual
orbitals to obtain orbitals better suited for correlated calculations.
The default option in DIRAC is to generate the modified virtual orbitals
by construction of a Fock operator for the system with some electrons
removed, as first proposed by Bauschlicher
\ :sup:``[1] <#_note-bauschlicher:mvo>`_`\ . The orbital string
specified with `.MVOVEC <#.MVOVEC>`_ indicates what occupied orbitals to
remove in the construction of the Hartree-Fock potential.

*Note* that the resulting modified virtual orbitals are no longer
canonical and can *not* be used in MP2 calculations; they can be used in
the `DIRRCI </dirac/index.php/Manual:DIRRCI>`_ and
`RELCCSD </dirac/index.php/Manual:RELCCSD>`_ modules.

.PRINT
------

Print level.

*Default:*

::

    .PRINT
     0

.VECMVO
-------

Read `orbital string </dirac/index.php/Manual:Orbital_strings>`_ of
orbitals to be included in the Hartree-Fock potential when constructing
the cationic Fock operator.

.INTFLG
-------

Specify what two-electron integrals to include during the construction
of the modified virtual orbitals (see
`.INTFLG </dirac/index.php/Manual:HAMILTONIAN#.INTFLG>`_ under
`\*\*HAMILTONIAN </dirac/index.php/Manual:HAMILTONIAN>`_).

*Default:* `.INTFLG </dirac/index.php/Manual:HAMILTONIAN#.INTFLG>`_ from
`\*\*HAMILTONIAN </dirac/index.php/Manual:HAMILTONIAN>`_.

References
----------

#. `↑ <#_ref-bauschlicher:mvo_0>`_ C. W. Bauschlicher, Jr., *The
   construction of modified virtual orbitals (MVO's) which are suited
   for configuration interaction calculations*, J. Chem. Phys. **72**,
   *880* (1980) `electronic
   version <http://dx.doi.org/10.1063/1.439243>`_.

