Manual:RELCCSD
==============

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

*Relativistic Coupled Cluster module RELCCSD* written by Lucas Visscher

This section specifies the input necessary for relativistic Coupled
Cluster calculations. The program is currently capable of doing energy
calculations at the MP2, CCSD, and CCSD(T) levels of theory
\ :sup:``[1] <#_note-visscher:relccsd1>`_`\ 
\ :sup:``[2] <#_note-pernpointner:relccsd3>`_`\ , while first order
properties are available only at the MP2 level of theory.

Open shell states can be calculated also via the Fock space Coupled
Cluster method \ :sup:``[3] <#_note-visscher:relccsd2>`_`\  in which
electrons are added to or subtracted from a closed shell core. All
sectors that involve no more than 2 replacements are available.

The input should be given in namelist form and is subdivided in a
general menu section and section describing energy, first-order,
second-order property and Fock-space Coupled Cluster calculations.

Contents
--------

-  `1 &RELCCSD <#.26RELCCSD>`_
-  `2 &CCENER <#.26CCENER>`_
-  `3 &CCSORT <#.26CCSORT>`_
-  `4 &CCFOPR <#.26CCFOPR>`_
-  `5 &CCSOPR <#.26CCSOPR>`_
-  `6 &CCFSPC <#.26CCFSPC>`_
-  `7 Sample inputs <#Sample_inputs>`_
-  `8 New input <#New_input>`_
-  `9 References <#References>`_

`&RELCCSD </dirac/index.php/Manual:RELCCSD:RELCCSD>`_
-----------------------------------------------------

Generic and control input.

`&CCENER </dirac/index.php/Manual:RELCCSD:CCENER>`_
---------------------------------------------------

Energy module.

`&CCSORT </dirac/index.php/Manual:RELCCSD:CCSORT>`_
---------------------------------------------------

Sorting module.

`&CCFOPR </dirac/index.php?title=Manual:RELCCSD:CCSFOPR&action=edit>`_
----------------------------------------------------------------------

First-order properties.

`&CCSOPR </dirac/index.php/Manual:RELCCSD:CCSOPR>`_
---------------------------------------------------

Second-order properties

`&CCFSPC </dirac/index.php/Manual:RELCCSD:CCFSPC>`_
---------------------------------------------------

Fock space Coupled Cluster module.

`Sample inputs </dirac/index.php/Manual:RELCCSD:Sample_inputs>`_
----------------------------------------------------------------

Some examples of the namelist input used in DIRAC10.

`New input </dirac/index.php/Manual:RELCCSD:New_input>`_
--------------------------------------------------------

The old namelist style input (described above) used in DIRAC10 will be
replaced by an input that is more consistent with the rest of the input
in the next release. This input is briefly described for those of you
who have access to the developer version of DIRAC.

References
----------

#. `↑ <#_ref-visscher:relccsd1_0>`_ L. Visscher, T. J. Lee, and K. G.
   Dyall, *Formulation and implementation of a relativistic unrestricted
   coupled cluster method including noniterative connected triples*, J.
   Chem. Phys. **105**, *8769* (1996) `electronic
   version <http://dx.doi.org/10.1063/1.472655>`_.
#. `↑ <#_ref-pernpointner:relccsd3_0>`_ M. Pernpointner and L. Visscher,
   *Parallelization of four-component calculations. II. Symmetry- driven
   parallelization of the 4-spinor CCSD algorithm*, J. Comput. Chem.
   **24**, *754* (2003) `electronic
   version <http://dx.doi.org/10.1002/jcc.10215>`_.
#. `↑ <#_ref-visscher:relccsd2_0>`_ L. Visscher, E. Eliav, and U.
   Kaldor, *Formulation and implementation of the relativistic
   Fock-space coupled cluster method for molecules*, J. Chem. Phys.
   **115**, *9720* (2001) `electronic
   version <http://dx.doi.org/10.1063/1.1415746>`_.

