:orphan:

.. _FDENMR:

================================================================
Calculation of NMR shieldings in Frozen Density Embedding scheme
================================================================

In this tutorial we look at calculations of NMR shielding tensor with London atomic orbitals
within the FDE scheme. 


Basic theory of NMR-LAO with FDE
================================

In order to calculate the NMR properties with LAOs, one must take into account 
modifications of the property gradient due to the use of perturbation-dependent basis sets.
Here we only present some key points, for the detailed explanations we encourage the reader
to consults papers...
The FDE contribution to the property gradient due to the use of perturbation-dependent basis sets:

.. math::

    E_{[1],emb}^{[B]} = - \int v_{emb}^I \Omega_{ia}^{B,I} - \iint w_{emb}^{I,I} \Omega_{ia}^{I} \Omega_{jj}^{I} - \iint w_{emb}^{I,II} \Omega_{ia}^{I} \Omega_{jj}^{II}

where :math:`\Omega_{pq}^{B}` is the first-order perturbed overlap distribution, which in a basis of London orbitals consists of two terms

.. math::

    \Omega_{pq}^{B} = \frac{i}{2} (R_{MN} \times r) \Omega_{pq} + (T_{pt}^{B*} \Omega_{tq} + \Omega_{pt}T_{tq}^B)

"direct" LAO term (the first term) and "reorthonormalization" term (two last terms in brackets).

In this tutorial we will show how to include each of the contributions to the property gradient in the
calculations of the NMR shielding tensor. We will also present how to calculate and plot the first-order
perturbed current density induced by magnetic field and how to use it in analysis of the NMR shielding tensor.

A sample calculation of NMR shielding tensor within FDE
=======================================================

We will consider a simple system of SeH\ :sub:`2`\ - H\ :sub:`2`\ O, in which SeH\ :sub:`2`\ 
constitutes an active subsystem (subsystem I) and H\ :sub:`2`\ O\ is an environment (subsystem II). 
The computational protocol is different depending on whether the coupling between two
subsystems is included or not (through the term dependent on :math:`w_{emb}^{I,II}`).
To avoid confusion, we will describe these two protocols separately.

We will need two molecular input files. File ``h2se.mol``, corresponding to an active subsystem I, reads:

.. literalinclude:: h2se.mol

and file ``h2o.mol``, standing for an environment subsystem II:

.. literalinclude:: h2o.mol

Calculations without the FDE coupling terms
-------------------------------------------

We will show how to include only the first two terms in the property gradient.

Step 1: SCF calculation on the subsystem II
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The first step is the calculation on subsystem II, the aim of which is
to get its density and electrostatic potential written on a text file.

The menu file ``getfrozen.inp`` reads:

.. literalinclude:: getfrozen.inp

and we make sure to recover the text file related to the :ref:`GRIDOU` keyword, which we name ``GRIDOUT``.
We should also remember to supply two files, 
related to keywords :ref:`GRIDOU` (``GRIDOUT``) and :ref:`EXONLY` (``FILEEX``),
which are copies of a numerical grid file

.. literalinclude:: getfrozen.run

After this calculation, the ``GRIDOUT`` file will be updated: additional columns on that file
correspond to the electrostatic potential and the density of subsystem II
calculated in each of grid points on that file.

Often, in FDE calculations, one uses the numerical grid which was generated for the supermolecular system.
If the numerical grid was generated in the DIRAC software, then it will be divided into batches of points.
However, FDE module is able to read only the continuous list of grid points, therefore such numerical grid has
to be reformatted. This can be done by the script, which can be found in ``utils`` directory of DIRAC:
``dft_to_fde_grid_convert.py``, which should be used on a ``numerical_grid`` file exported from DIRAC:

.. literalinclude:: grid_convert

The original grid file is saved as ``numerical_grid.ORIG``

Step 2: SCF and response calculation on the subsystem I
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Now we can use the information on the subsystem II written on the ``GRIDOUT`` file for constructing
the embedding potential, :math:`v_{emb}^{I}`, and embedding kernel, :math:`w_{emb}^{I,I}`.
In order to do that, we use the ``GRIDOUT`` file saved from Step 1 and copy it to a ``FROZEN`` file,
given as an argument to the :ref:`FRDENS` keyword.

To simplify this tutorial, we will do the SCF step and response calculations in one run, by using the following input file:

.. literalinclude:: updatefrozen_v.inp

and the ``pam`` command is:

.. literalinclude:: updatefrozen_v.run

As we used the keyword :ref:`LAODPT`, the property gradient will contain the contribution
from the embedding potential, :math:`v_{emb}^{I}`.
In order to include also contributions from the embedding kernel, :math:`w_{emb}^{I,I}`, 
we should use two additional keywords in the :ref:`*FDE` Section:

.. literalinclude:: updatefrozen_v_w11.inp

with the ``pam`` command:

.. literalinclude:: updatefrozen_v_w11.inp

Calculations with the FDE coupling terms
----------------------------------------

We will show how to include the coupling between two subsystems in the property gradient 
(term dependent on :math:`w_{emb}^{I,II}`).

Step 1: SCF and response calculation on the subsystem II
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The first step is the calculation on subsystem II, the aim of which is
to get its density and electrostatic potential written on a text file.
In addition, we will export the perturbed density of subsystem II

The menu file ``getfrozenpert.inp`` reads:

.. literalinclude:: getfrozenpert.inp


...TODO gosia


A sample calculation of magnetically-induced current density within FDE
=======================================================================

Here we will demonstrate how to calculate and visualize the first-order current density induced by magnetic field (see :cite:`Bast2009vis`):

.. math::

    j^{B} = -ec \sum_i \left\(\psi_i^B \alpha \psi_i  + \psi_i \alpha \psi_i^B\right)

with the first-order perturbed molecular orbitals defined as (see :cite:`Ilias2009`):

.. math::

    \psi_i^{B} = c_{\mu i}\omega_{\mu}^B + \psi_p T_{pi}^B -\psi_{a}\kappa_{ai}^B

We will also show how to visualize the NMR shielding density:

.. math::

    \sigma^K = -\frac{1}{c^2} \int dr \frac{1}{r_K^3} (r_K \times j^B)

with the use of the magnetically-induced current density.
We will focus on these properties for subsystem I.

In order to get an access to these properties, we should perform SCF and response calculations,
in steps as described above. The only change is to retrieve appropriate files from the final response 
calculations on subsystem I: ``PAMXVC``, ``TBMO`` and ``DFCOEF`` files.
The ``pam`` command should be following:

.. literalinclude:: response_vis.run

First, we should identify the response vectors on the ``PAMXVC`` file, by using the following command:

.. literalinclude:: pamxvc_cmd

where ``rspread.x`` can be found in the DIRAC's ``build`` directory.
An example output could be analogous to this one:

..TODO

The input which allows to calculate and plot the magnetically-induced current density can be as follows:

.. literalinclude:: vis_j.inp

where the perturbed current density will be calculated in a plane of SeH2 molecule (for full list of possibilities see :ref:`**VISUAL` and the VISUAL - tutorial).
And the corresponding ``pam`` command:

.. literalinclude:: vis_j.run

From this calculation, we retrieve files ``plot.*``, which can be visualized in any plotting software.
In order to see the effect of the use of FDE, one may substract the plot calculated for SeH2 with environment accounted
for by FDE scheme from the plot calculated for an isolated SeH2 molecule.
As a result we get:

..TODO






