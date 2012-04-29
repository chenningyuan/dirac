Manual:PROPERTIES
=================

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

Contents
--------

-  `1 General control statements <#General_control_statements>`_

   -  `1.1 .PRINT <#.PRINT>`_
   -  `1.2 .ABUNDANCIES <#.ABUNDANCIES>`_
   -  `1.3 .RKBIMP <#.RKBIMP>`_

-  `2 Predefined electric properties <#Predefined_electric_properties>`_

   -  `2.1 .DIPOLE <#.DIPOLE>`_
   -  `2.2 .QUADRUPOLE <#.QUADRUPOLE>`_
   -  `2.3 .EFG <#.EFG>`_
   -  `2.4 .NQCC <#.NQCC>`_
   -  `2.5 .POLARIZABILITY <#.POLARIZABILITY>`_
   -  `2.6 .FIRST ORDER
      HYPERPOLARIZABILITY <#.FIRST_ORDER_HYPERPOLARIZABILITY>`_
   -  `2.7 .VERDET <#.VERDET>`_
   -  `2.8 .TWO-PHOTON <#.TWO-PHOTON>`_

-  `3 Predefined magnetic properties <#Predefined_magnetic_properties>`_

   -  `3.1 .NMR <#.NMR>`_
   -  `3.2 .SHIELDING <#.SHIELDING>`_
   -  `3.3 .MAGNET <#.MAGNET>`_
   -  `3.4 .SPIN-SPIN COUPLING <#.SPIN-SPIN_COUPLING>`_
   -  `3.5 .DSO <#.DSO>`_
   -  `3.6 .NSTDIAMAGNETIC <#.NSTDIAMAGNETIC>`_

-  `4 Other predefined properties <#Other_predefined_properties>`_

   -  `4.1 .MOLGRD <#.MOLGRD>`_
   -  `4.2 .PVC <#.PVC>`_
   -  `4.3 .RHONUC <#.RHONUC>`_
   -  `4.4 .EFFDEN <#.EFFDEN>`_

-  `5 References <#References>`_

This section allows for the evaluation of a large number of molecular
properties. Available properties include:

-  Expectation values (e.g. dipole moment and electric field gradients).
-  Linear response properties (e.g. polarizability and NMR parameters).
-  Quadratic response properties (e.g. hyperpolarizabilities).

For convenience some common properties can be specified directly in this
section, which means that the user in principle does not need to know
how they are calculated. Note, however, that response functions are by
default static, but frequencies can be added in the relevant
subsections.

Properties which are not predefined must be specified in detail in the
relevant input section (see the `specification of one-electron
operators </dirac/index.php/Manual:One-electron_operators>`_).

By default no properties are calculated.

General control statements
==========================

.PRINT
------

Print level.

*Default:*

::

    .PRINT
     0

.ABUNDANCIES
------------

For properties that make reference to isotopes, give threshold level (in
% abundance) for isotopes to print.

*Default:*

::

    .ABUNDANCIES
     1.0

.RKBIMP
-------

Import coefficients calculated with restricted kinetic balance (RKB) in
a calculation using unrestricted kinetic balance (UKB). This option is a
simple way to generated restricted magnetic balance for the calculation
of NMR shieldings. This option works in the general SO case, but not in
the spinfree case since spinfree calculations are not possible with UKB.

Predefined electric properties
==============================

.DIPOLE
-------

Evaluate the electronic dipole moment (expectation value).

.QUADRUPOLE
-----------

Evaluate the electronic quadrupole moment (expectation value).

.EFG
----

Evaluate electric field gradients (expectation values).

Atomic centers may be restricted with
`.SELECT </dirac/index.php/Manual:INTEGRALS#.SELECT>`_ under
`\*\*INTEGRALS </dirac/index.php/Manual:INTEGRALS>`_.

.NQCC
-----

Evaluate nuclear quadrupole coupling constants (expectation values).

Atomic centers may be restricted with
`.SELECT </dirac/index.php/Manual:INTEGRALS#.SELECT>`_ under
`\*\*INTEGRALS </dirac/index.php/Manual:INTEGRALS>`_.

.POLARIZABILITY
---------------

Evaluate the static electronic dipole polarizability tensor (linear
response).

.FIRST ORDER HYPERPOLARIZABILITY
--------------------------------

Evaluate static electronic dipole first-order hyperpolarizability tensor
(quadratic response). Results are also given for the static electronic
dipole polarizability.

.VERDET
-------

Evaluate Verdet constants\ :sup:``[1] <#_note-ekstrom:verdet>`_`\ 
(quadratic response) for a dynamic electric field corresponding to Ruby
laser wavelength of 694 nm and a static magnetic field along the
propagation direction of the light beam (in this case, the default
frequencies of the quadratic response function thus become
ω\ :sub:`*B*`\  = 0.0656 and ω\ :sub:`*C*`\  = 0.0). A Verdet
calculation cannot be specified in combination with other quadratic
response calculations.

The frequencies can be changed using
`.BFREQ </dirac/index.php/Manual:PROPERTIES:QUADRATIC_RESPONSE#.BFREQ>`_
in `\*QUADRATIC
RESPONSE </dirac/index.php/Manual:PROPERTIES:QUADRATIC_RESPONSE>`_.

.TWO-PHOTON
-----------

Evaluate two-photon absorption cross sections (quadratic response). Give
the number of desired states in each boson symmetry. Cannot be specified
in combination with other quadratic response calculations.

*Example:* Point group with four boson irreps, (e.g.
*C*\ :sub:`2*v*`\ ).

::

    .TWO-PHOTON
     5 5 5 0

Predefined magnetic properties
==============================

.NMR
----

Evaluate nuclear magnetic shieldings and indirect spin-spin couplings
(linear response).

Atomic centers may be restricted with
`.SELECT </dirac/index.php/Manual:INTEGRALS#.SELECT>`_ under
`\*\*INTEGRALS </dirac/index.php/Manual:INTEGRALS>`_.

.SHIELDING
----------

Evaluate nuclear magnetic shieldings (linear response). `Print
level <#.PRINT>`_ 2 gives tensor and longer output. `Print
level <#.PRINT>`_ 4 gives the raw values in symmetry coordinates as
well.

Atomic centers may be restricted with
`.SELECT </dirac/index.php/Manual:INTEGRALS#.SELECT>`_ under
`\*\*INTEGRALS </dirac/index.php/Manual:INTEGRALS>`_.

.MAGNET
-------

Evaluate the magnetic susceptibilities tensor (linear response and
expectation values).

.SPIN-SPIN COUPLING
-------------------

Evaluate indirect spin-spin couplings (linear response).

Atomic centers may be restricted with
`.SELECT </dirac/index.php/Manual:INTEGRALS#.SELECT>`_ under
`\*\*INTEGRALS </dirac/index.php/Manual:INTEGRALS>`_.

.DSO
----

Evaluate the diamagnetic contribution to indirect spin-spin couplings as
an expectation value.

Atomic centers may be restricted with
`.SELECT </dirac/index.php/Manual:INTEGRALS#.SELECT>`_ under
`\*\*INTEGRALS </dirac/index.php/Manual:INTEGRALS>`_.

.NSTDIAMAGNETIC
---------------

Evaluate the diamagnetic contribution to nuclear magnetic shieldings as
an expectation value.

Atomic centers may be restricted with
`.SELECT </dirac/index.php/Manual:INTEGRALS#.SELECT>`_ under
`\*\*INTEGRALS </dirac/index.php/Manual:INTEGRALS>`_.

Other predefined properties
===========================

.MOLGRD
-------

Evaluate the molecular gradient, i.e. |{\\frac{{\\partial E}}{{\\partial
X\_A}}}| where *X*\ :sub:`*A*`\  are the coordinates of the nuclei. This
is an expectation value of one- and two-electron operators. Normally the
molecular gradient evaluation is not invoked explicitly with this
keyword but rather implicitly in the geometry optimization module.

.PVC
----

Calculate matrix elements over the nuclear spin-independent
parity-violating operator, e.g. calculate energy differences between
enantiomers \ :sup:``[2] <#_note-laerdahl:epv>`_`\ .

.RHONUC
-------

Calculate electronic density at the nucleus(nuclei).

Atomic centers may be restricted with
`.SELECT </dirac/index.php/Manual:INTEGRALS#.SELECT>`_ under
`\*\*INTEGRALS </dirac/index.php/Manual:INTEGRALS>`_.

.EFFDEN
-------

Calculate effective electronic density at the
nucleus(nuclei).\ :sup:``[3] <#_note-sknecht:mossbauer>`_`\ .

Atomic centers may be restricted with
`.SELECT </dirac/index.php/Manual:INTEGRALS#.SELECT>`_ under
`\*\*INTEGRALS </dirac/index.php/Manual:INTEGRALS>`_.

References
==========

#. `↑ <#_ref-ekstrom:verdet_0>`_ Ulf Ekstrom and Patrick Norman and
   Antonio Rizzo, *Four-component Hartree--Fock calculations of
   magnetic-field induced circular birefringence---Faraday effect---in
   noble gases and dihalogens*, The Journal of Chemical Physics **122**,
   *074321* (2005) `electronic
   version <http://dx.doi.org/10.1063/1.1849167>`_.
#. `↑ <#_ref-laerdahl:epv_0>`_ J. K. Lærdahl and P. Schwerdtfeger,
   *Fully relativistic ab initio calculations of the energies of chiral
   molecules including partiy-violating weak interactions*, Phys. Rev. A
   **60**, *4439* (1999) `electronic
   version <http://dx.doi.org/10.1103/PhysRevA.60.4439>`_.
#. `↑ <#_ref-sknecht:mossbauer_0>`_ S. Knecht, S. Fux, R. van Meer, L.
   Visscher, M. Reiher, T. Saue, *Mössbauer spectroscopy for heavy
   elements: a relativistic benchmark study of mercury*, Theor. Chem.
   Acc. **129**, *631-650* (2011) `electronic
   version <http://dx.doi.org/10.1007/s00214-011-0911-2>`_.

Retrieved from
"`http://wiki.chem.vu.nl/dirac/index.php/Manual:PROPERTIES <http://wiki.chem.vu.nl/dirac/index.php/Manual:PROPERTIES>`_"

Views
'''''

-  `Article </dirac/index.php/Manual:PROPERTIES>`_
-  `Discussion </dirac/index.php?title=Talk:Manual:PROPERTIES&action=edit>`_
-  `View source </dirac/index.php?title=Manual:PROPERTIES&action=edit>`_
-  `History </dirac/index.php?title=Manual:PROPERTIES&action=history>`_

Personal tools
''''''''''''''

-  `Log in / create
   account </dirac/index.php?title=Special:Userlogin&returnto=Manual:PROPERTIES>`_

` </dirac/index.php/Dirac_Program>`_

Navigation
''''''''''

-  `Home </dirac/index.php/Dirac_Program>`_
-  `Documentation </dirac/index.php/Documentation>`_
-  `Reference manual </dirac/index.php/Manual>`_
-  `Discussion forum </dirac/index.php/DiscussionForum>`_
-  `News Archive </dirac/index.php/News_archive>`_
-  `Basis sets </dirac/index.php/Basis_Sets>`_

Search
''''''

Toolbox
'''''''

-  `What links
   here </dirac/index.php/Special:Whatlinkshere/Manual:PROPERTIES>`_
-  `Related
   changes </dirac/index.php/Special:Recentchangeslinked/Manual:PROPERTIES>`_
-  `Upload file </dirac/index.php/Special:Upload>`_
-  `Special pages </dirac/index.php/Special:Specialpages>`_
-  `Printable
   version </dirac/index.php?title=Manual:PROPERTIES&printable=yes>`_
-  `Permanent
   link </dirac/index.php?title=Manual:PROPERTIES&oldid=30141>`_

`|Powered by MediaWiki| <http://www.mediawiki.org/>`_

-  This page was last modified 12:15, 2 August 2011.
-  This page has been accessed 2,450 times.
-  `Privacy policy </dirac/index.php/Diracwiki:Privacy_policy>`_
-  `About Diracwiki </dirac/index.php/Diracwiki:About>`_
-  `Disclaimers </dirac/index.php/Diracwiki:General_disclaimer>`_

.. |{\\frac{{\\partial E}}{{\\partial
X\_A}}}| image:: /dirac/images/math/2/0/c/20c2e6537efd4db350d7545b3b1147f7.png
.. |Powered by
MediaWiki| image:: /dirac/skins/common/images/poweredby_mediawiki_88x31.png
