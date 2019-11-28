:orphan:
 

starstar(PROPERTIES)

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

keyword(PRINT)

Print level.

*Default:*

::

    .PRINT
     0

keyword(ABUNDANCIES)

For properties that make reference to isotopes, give threshold level (in
% abundance) for isotopes to print.

*Default:*

::

    .ABUNDANCIES
     1.0

keyword(RKBIMP)

Import coefficients calculated with restricted kinetic balance (RKB) in
a calculation using unrestricted kinetic balance (UKB). This option is a
simple way to generated restricted magnetic balance for the calculation
of NMR shieldings. This option works in the general SO case, but not in
the spinfree case since spinfree calculations are not possible with UKB.

keyword(NOPCTR)

In two-component infinite-order relativistic calculations (with :ref:`HAMILTONIAN_.X2C`) take only LL block
of four-component property operators to avoid the picture change transformation.
Experimental option, use with care.


Predefined electric properties
==============================

keyword(DIPOLE)

Evaluate the electronic dipole moment (expectation value).

keyword(QUADRUPOLE)

Evaluate the electronic quadrupole moment (expectation value).

keyword(EFG)

Evaluate electric field gradients (expectation values).

Atomic centers may be restricted with :ref:`INTEGRALS_.SELECT` under :ref:`**INTEGRALS`.

keyword(NQCC)

Evaluate nuclear quadrupole coupling constants (expectation values).

Atomic centers may be restricted with :ref:`INTEGRALS_.SELECT` under :ref:`**INTEGRALS`.

keyword(POLARIZABILITY)

Evaluate the static electronic dipole polarizability tensor (linear
response).

keyword(FIRST ORDER HYPERPOLARIZABILITY)

Evaluate static electronic dipole first-order hyperpolarizability tensor
(quadratic response). Results are also given for the static electronic
dipole polarizability.

keyword(VERDET)

Evaluate Verdet constants :cite:`Ekstrom2005`
(quadratic response) for a dynamic electric field corresponding to Ruby
laser wavelength of 694 nm and a static magnetic field along the
propagation direction of the light beam (in this case, the default
frequencies of the quadratic response function thus become
ω\ :sub:`*B*`\  = 0.0656 and ω\ :sub:`*C*`\  = 0.0). A Verdet
calculation cannot be specified in combination with other quadratic
response calculations.

The frequencies can be changed using :ref:`QUADRATIC_RESPONSE_.B FREQ` in :ref:`*QUADRATIC RESPONSE`.

keyword(TWO-PHOTON)

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

keyword(NMR)

Evaluate nuclear magnetic shieldings and indirect spin-spin couplings
(linear response).

Atomic centers may be restricted with :ref:`INTEGRALS_.SELECT` under :ref:`**INTEGRALS`.

See below for advice on calculation of diamagnetic terms.

keyword(SHIELDING)

Evaluate nuclear magnetic shieldings (linear response). :ref:`PROPERTIES_.PRINT`
2 gives tensor and longer output. The :ref:`PROPERTIES_.PRINT` 4
gives the raw values in symmetry coordinates as well.

Atomic centers may be restricted with :ref:`INTEGRALS_.SELECT` under :ref:`**INTEGRALS`.

keyword(MAGNET)

Evaluate the magnetic susceptibilities tensor (linear response and
expectation values).

keyword(ROTG)

Evaluate rotational g-tensors: linear response and nuclear contributions.

.. note::

   The current implementation gives theoretically correct results only for molecules in equilibrium.

   
Results are dimensionless.

The total g-tensor, as well as its linear response and nuclear contributions are always given separately.

Using :ref:`PROPERTIES_.PRINT` 2, the paramagnetic (e-e) and diamagnetic (e-p) parts of the linear response contributions are given separately, together with results for the :math:`\mathbf{L}` and :math:`\mathbf{S}` parts of the linear response.

See :cite:`Aucar_JCP2014` for further details.

keyword(SPIN-SPIN COUPLING)

Evaluate indirect spin-spin couplings (linear response).

Atomic centers may be restricted with :ref:`INTEGRALS_.SELECT` under :ref:`**INTEGRALS`.

The default is to calculate the diamagnetic term via occupied positive energy to virtual negative energy orbital rotations
(also called electron-positron rotations), see :cite:`Aucar1999` for the theory.
The quality of this is very basis set dependent.
It is generally more accurate to use the non-relativistic expectation value expression
for the diamagnetic term, activated with keyword .DSO in this section.
You must also add :ref:`LINEAR_RESPONSE_.SKIPEP` under :ref:`*LINEAR RESPONSE` to 
exclude the diamagnetic term from the linear response calculation.

keyword(DSO)

Evaluate the diamagnetic contribution to indirect spin-spin couplings as
an expectation value of the non-relativistic DSO operator.

Atomic centers may be restricted with :ref:`INTEGRALS_.SELECT` under :ref:`**INTEGRALS`.

keyword(NSTDIAMAGNETIC)

Evaluate the diamagnetic contribution to nuclear magnetic shielding tensor as
an expectation value of the non-relativistic operator.

Atomic centers may be restricted with :ref:`INTEGRALS_.SELECT` under :ref:`**INTEGRALS`.

Other predefined properties
===========================

keyword(MOLGRD)

Evaluate the molecular gradient, i.e.

.. math::

  \frac{\partial E}{\partial \mathbf{X}_A}

where :math:`\mathbf{X}_{A}` are the coordinates of the nuclei. This
is an expectation value of one- and two-electron operators. Normally the
molecular gradient evaluation is not invoked explicitly with this
keyword but rather implicitly in the geometry optimization module.

keyword(PVC)

Calculate matrix elements over the nuclear spin-independent
parity-violating operator, e.g. calculate energy differences between
enantiomers, see :cite:`Laerdahl1999` .

keyword(RHONUC)

Calculate electronic density at the nucleus(nuclei).

Atomic centers may be restricted with :ref:`INTEGRALS_.SELECT` under :ref:`**INTEGRALS`.

keyword(EFFDEN)

Calculate effective electronic density at the
nucleus(nuclei), :cite:`Knecht2011`.

Atomic centers may be restricted with :ref:`INTEGRALS_.SELECT` under :ref:`**INTEGRALS`.

keyword(SPIN-ROTATION)

Evaluate nuclear spin-rotation constants: linear response, expectation value and nuclear contributions.

.. note::

   The current implementation gives theoretically correct results only for molecules in equilibrium.

   
Results are given in kHz, but for particular :ref:`PROPERTIES_.PRINT` values they could also be given in ppm (to compare results with :ref:`PROPERTIES_.SHIELDING`).

The total spin-rotation tensors, as well as their linear response, expectation value and nuclear contributions are always given separately.

Using :ref:`PROPERTIES_.PRINT` 0 or 1, results are only given in kHz, whereas employing :ref:`PROPERTIES_.PRINT` 2 or 3 they are also shown in ppm.

In addition, when :ref:`PROPERTIES_.PRINT` 1 or 3 are used, the paramagnetic-like (e-e) and diamagnetic-like (e-p) parts of the linear response contributions are given separately, together with results for the :math:`\mathbf{L}` and :math:`\mathbf{S}` parts of the linear response.

Finally, employing :ref:`PROPERTIES_.PRINT` 4, nuclear and expectation value contributions are given separately (with Thomas precesion included). More details in the Tutorial Section.

See :cite:`Aucar_JCP2012` for further details.

