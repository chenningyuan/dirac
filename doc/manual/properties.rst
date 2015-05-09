
.. _properties:

==============
\*\*PROPERTIES
==============

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
relevant input section (see :ref:`one_electron_operators`).

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

Atomic centers may be restricted with :ref:`.SELECT` under :ref:`**INTEGRALS`.

.NQCC
-----

Evaluate nuclear quadrupole coupling constants (expectation values).

Atomic centers may be restricted with :ref:`.SELECT` under :ref:`**INTEGRALS`.

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

Evaluate Verdet constants :cite:`Ekstrom2005`
(quadratic response) for a dynamic electric field corresponding to Ruby
laser wavelength of 694 nm and a static magnetic field along the
propagation direction of the light beam (in this case, the default
frequencies of the quadratic response function thus become
ω\ :sub:`*B*`\  = 0.0656 and ω\ :sub:`*C*`\  = 0.0). A Verdet
calculation cannot be specified in combination with other quadratic
response calculations.

The frequencies can be changed using :ref:`.BFREQ` in :ref:`*QUADRATIC RESPONSE`.

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

Atomic centers may be restricted with :ref:`.SELECT` under :ref:`**INTEGRALS`.

.SHIELDING
----------

Evaluate nuclear magnetic shieldings (linear response). :ref:`.PRINT`
2 gives tensor and longer output. :ref:`.PRINT` 4
gives the raw values in symmetry coordinates as
well.

Atomic centers may be restricted with :ref:`.SELECT` under :ref:`**INTEGRALS`.

.MAGNET
-------

Evaluate the magnetic susceptibilities tensor (linear response and
expectation values).

.SPIN-SPIN COUPLING
-------------------

Evaluate indirect spin-spin couplings (linear response).

Atomic centers may be restricted with :ref:`.SELECT` under :ref:`**INTEGRALS`.

.DSO
----

Evaluate the diamagnetic contribution to indirect spin-spin couplings as
an expectation value.

Atomic centers may be restricted with :ref:`.SELECT` under :ref:`**INTEGRALS`.

.NSTDIAMAGNETIC
---------------

Evaluate the diamagnetic contribution to nuclear magnetic shieldings as
an expectation value.

Atomic centers may be restricted with :ref:`.SELECT` under :ref:`**INTEGRALS`.

Other predefined properties
===========================

.MOLGRD
-------

Evaluate the molecular gradient, i.e.

..|{\\frac{{\\partial E}}{{\\partial X\_A}}}|

where *X*\ :sub:`*A*`\  are the coordinates of the nuclei. This
is an expectation value of one- and two-electron operators. Normally the
molecular gradient evaluation is not invoked explicitly with this
keyword but rather implicitly in the geometry optimization module.

.PVC
----

Calculate matrix elements over the nuclear spin-independent
parity-violating operator, e.g. calculate energy differences between
enantiomers, see :cite:`Laerdahl1999` .

.RHONUC
-------

Calculate electronic density at the nucleus(nuclei).

Atomic centers may be restricted with :ref:`.SELECT` under :ref:`**INTEGRALS`.

.EFFDEN
-------

Calculate effective electronic density at the
nucleus(nuclei), :cite:`Knecht2011`.

Atomic centers may be restricted with :ref:`.SELECT` under :ref:`**INTEGRALS`.



