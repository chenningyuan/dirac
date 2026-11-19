:orphan:

========================================================
Parity-violation contribution to electric field gradient
========================================================

Introduction
------------

In this tutorial we introduce the calculation of parity violation contribution to the electric field gradient
tensor as given in the DIRAC code.

As described in :doc:`PROPERTIES </manual/properties>`,
the parity-violation contribution to the electronic part of the electric field gradient tensor (PVCEFG) of a nucleus :math:`K` is given by

.. math::

   q^{PV}_{ij} (\mathbf{R}_K) =  \langle\langle \hat q_{ij} (\mathbf{R}_K);     \hat{H}^{PV}_K \rangle\rangle_0

where the contribution from the K nucleus to the nuclear spin-independent parity-violating operator

.. math::

   \hat{H}^{PV}_K =\frac{G_F}{2\sqrt{2}}\sum_{i,K}Q_{w,K}\gamma_i^5\rho_K(\mathbf{r}_i)

is expressed in terms of the Fermi coupling constant :math:`G_F=2.222516\times 10^{-14}E_ha_0^3` (this value corresponds to CODATA 2022, but it can change if you use other CODATA set),
the weak nuclear charge :math:`Q_{w,K}` and the normalized nuclear charge density :math:`\rho_K` (in units of the inverse of cube distances).
Results are given in a.u. Further details can be found in :cite:`Aucar2025_PVEFG`.


Application to the H2O2 molecule
--------------------------------

As an example, we show a calculation of the PV contribution to the EFG tensor of the oxygen and hydrogen nuclei at the Hydrogen peroxide molecule. The input file `efg_dc.inp` is given by

.. literalinclude:: efg_dc.inp

whereas the molecular input file `H2O2.mol` is

.. literalinclude:: H2O2.mol

The calculation is run using::

   pam --inp=efg_dc --mol=H2O2


As a result, PVCEFG are obtained at the RPA (i.e., coupled Hartree-Fock) level of approach. The code also works at the DFT level.

Reading the output file
-----------------------

As the :ref:`PROPERTIES_.PRINT` flag in the input file is set to 1, extra details are printed. The results are given in au.
The PVCEFG to the oxygen and hydrogen nuclei will look like:

.. literalinclude:: PVCEFG-H2O2

As it can be seen, the total PVC to the EFG tensor elements is further separated in their (e-e) and (e-p) parts. 


..  warning:: The basis set used here is rather small for demonstration purposes; larger basis sets are recommended for production calculations.