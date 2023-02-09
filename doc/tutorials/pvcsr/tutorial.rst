:orphan:

================================================================
Parity-violation contribution to nuclear spin-rotation constants
================================================================

Introduction
------------

In this tutorial we introduce the calculation of parity violation contribution to nuclear spin-rotation tensors
as given in the DIRAC code.

The parity-violation contribution to nuclear spin-rotation (PVCNSR) tensor of a nucleus :math:`N` is given by

.. math::

   {\bf M}_N^{PV} = \; \frac{-\hslash}{h} \frac{G_F \, (1-4sen^2\theta_W)}{\sqrt{2} \, c} \, \lambda_N  \;  \langle \langle \rho_N({\bf r}) \, c \, {\bf \alpha} \; ; \; {\bf J}_e\rangle\rangle \; \cdot  \; {\bf I}^{-1}

where :math:`\rho_N({\bf r})` is the normalized nuclear charge density; :math:`{\bf I}` is the inertia tensor of the molecule and :math:`{\bf J}_e = \left({\bf r} - {\bf R}_{CM} \right) \times {\bf p}+{\bf S}_e` is the electronic total angular momentum.


Application to the H2O2 molecule
--------------------------------

As an example, we show a calculation of the PV contribution to the NSR constant of the oxygen nucleus at the Hydrogen peroxide molecule. The input file `pvcsr.inp` is given by

.. literalinclude:: pvcsr.inp

whereas the molecular input file `H2O2.mol` is

.. literalinclude:: H2O2.mol

The calculation is run using::

   pam --inp=pvcsr --mol=H2O2


As a result, PVCNSR are obtained at the coupled Hartree-Fock level. The code also works at the DFT level.


Reading the output file
-----------------------

As the :ref:`PROPERTIES_.PRINT` flag in the input file is set to 1, the results are fully detailed and given in Hz.

The PVCNSR of the oxygen nucleus will look like:

.. literalinclude:: PVCNSR-H2O2

As it can be seen, the total SR constant of the oxygen nucleus is given. The linear response function is further separated in their (e-e) and (e-p) parts, as well as in their :math:`\mathbf{L}` and :math:`\mathbf{S}` parts.
