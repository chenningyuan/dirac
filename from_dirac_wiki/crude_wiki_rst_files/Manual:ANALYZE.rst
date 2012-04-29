Manual:ANALYZE
==============

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

\*\*ANALYZE -- analyzing tools
==============================

This section is aimed to analyze the final Hartree-Fock wave function by
using one or more analysis modules. By default none of them is
activated.

Contents
--------

-  `1 \*\*ANALYZE -- analyzing
   tools <#.2A.2AANALYZE_--_analyzing_tools>`_

   -  `1.1 .PRIVEC <#.PRIVEC>`_
   -  `1.2 .MULPOP <#.MULPOP>`_
   -  `1.3 .PROJEC <#.PROJEC>`_
   -  `1.4 .DENSIT <#.DENSIT>`_
   -  `1.5 .LOCALI <#.LOCALI>`_
   -  `1.6 References <#References>`_

.PRIVEC
-------

Print vectors. Activates the
`\*PRIVEC </dirac/index.php/Manual:ANALYZE:PRIVEC>`_ subsection.

.MULPOP
-------

Perform Mulliken population analysis.
\ :sup:``[1] <#_note-mulliken-popana>`_`\ . Continues to
`\*MULPOP </dirac/index.php/Manual:ANALYZE:MULPOP>`_ subsection.

.PROJEC
-------

Perform projection analysis. \ :sup:``[2] <#_note-saue-faegri-01>`_`\ .
Activates the
`\*PROJECTION </dirac/index.php/Manual:ANALYZE:PROJECTION>`_ subsection.

.DENSIT
-------

Write density to a formatted file in Gaussian cube format. Activates the
`\*DENSITY </dirac/index.php/Manual:ANALYZE:DENSITY>`_ subsection.

.LOCALI
-------

Localize orbitals using the Pipek-Mezey criterion.
\ :sup:``[3] <#_note-dubillard-saue-06>`_`\ . Continues to the activated
`\*LOCALIZATION </dirac/index.php/Manual:ANALYZE:LOCALIZATION>`_
subsection. Note that in the present implementation this only works in
C\ :sub:`1`\  symmetry.

References
----------

#. `↑ <#_ref-mulliken-popana_0>`_ R. S. Mulliken, *Electronic Population
   Analysis on LCAO-MO Molecular Wave Functions. I*, J. Chem. Phys.
   **23**, *1833* (1955) `electronic
   version <http://dx.doi.org/10.1063/1.1740588>`_.
#. `↑ <#_ref-saue-faegri-01_0>`_ K.Fægri and T. Saue, *Diatomic
   molecules between very heavy elements of group 13 and group 17 - a
   study of relativistic effects on bonding*, J. Chem. Phys. **115**,
   *2456* (2001) `electronic
   version <http://dx.doi.org/10.1063/1.1385366>`_.
#. `↑ <#_ref-dubillard-saue-06_0>`_ S. Dubillard, J.-B. Rota, T. Saue
   and K.Fægri, *Bonding analysis using localized relativistic orbitals:
   Water, the ultrarelativistic case and the heavy homologues
   H\ :sub:`2`\ X (X=Te, Po, eka-Po)*, J. Chem. Phys. **124**, *154307*
   (2006) `electronic version <http://dx.doi.org/10.1063/1.2187001>`_.

