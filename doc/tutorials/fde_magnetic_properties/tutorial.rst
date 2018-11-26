:orphan:

.. _FDEMAG:


=============================================================================================================
NMR shielding and magnetizability in Frozen Density Embedding (FDE) scheme with London atomic orbitals (LAOs)
=============================================================================================================

In this tutorial we will show how to calculate the NMR shielding tensor and the magnetizability tensor
with London atomic orbitals within the FDE scheme. 



NMR shielding and magnetizability in FDE scheme - theory overview
=================================================================

The NMR shielding tensor and the magnetizability tensor (which we will jointly denote as :math:`P`) 
can be calculated as a sum of a linear response contribution (:math:`P_{linres}`) 
and an expectation value term (:math:`P_{expval}`):

.. math::
   :label: P

    P = P_{linres} + P_{expval}

The former contribution is determined from the response equations, which require the modification of the
property gradient due to the use of perturbation-dependent basis set in the FDE scheme (discussed in the next paragraph). 
The expectation value contribution requires analogous modifications only for the magnetizability tensor
and will be discussed in a separate paragraph.
 


Modifications of the property gradient
--------------------------------------

In order to calculate the NMR shielding tensor or the magnetizability tensor with LAOs, one must take into account 
modifications of the property gradient due to the use of perturbation-dependent basis sets (see :cite:`Ilias2009`).
In FDE scheme, further modifications are needed, since the property gradient includes not only terms from
isolated subsystems (here marked as "I" and "II"), 
but also a contribution from the interaction energy between these subsystems (see :cite:`Hofener2013`).

For instance, these new contributions to the property gradient of subsystem "I" can be written as (see :cite: `olejniczak2017`):

.. math::
   :label: fde_pg

    E_{[1],emb}^{[B],I} = - \int v_{emb}^I \Omega_{ia}^{B,I} - \iint w_{emb}^{I,I} \Omega_{ia}^{I} \Omega_{jj}^{B,I} - \iint w_{emb}^{I,II} \Omega_{ia}^{I} \Omega_{jj}^{B,II}

where :math:`\Omega_{pq}^{B,X}` is the first-order perturbed overlap distribution for a subsystem X (:math:`X\in\{I,II\}`), 
which in a basis of London orbitals consists of two terms

.. math::
   :label: mo_lao_b

    \Omega_{pq}^{B,X} = \frac{i}{2} (R_{MN} \times r) \Omega_{pq}^X + (T_{pt}^{B*,X} \Omega_{tq}^X + \Omega_{pt}^XT_{tq}^{B,X}),

"direct" LAO term (the first term) and "reorthonormalization" term (in parenthesis).
The :math:`v_{emb}^I` is the embedding potential of subsystem I, :math:`w_{emb}^{I,I}` is the embedding kernel involving subsystem I 
and :math:`w_{emb}^{I,II}` is the embedding kernel involving the coupling between subsystem I and subsystem II.

In this tutorial we will show how to include each of these contributions to the property gradient in the
calculations of the NMR shielding tensor or the magnetizability tensor. 



Modifications of the expectation value contributions to the magnetizability tensor
----------------------------------------------------------------------------------

While the :math:`E_{[1],emb}^{[B],X}` contribution to the property gradient formulated in response equations of subsystem "X" is the
only modification that needs to be accounted for in the calculations of NMR shieldings, the magnetizability tensor requires
additional terms in the expectation value part which arise due to the use of perturbation-dependent basis sets.
As the magnetizability tensor is an extensive property, its expectation value part reads

.. math::
   :label: fde_magn_expval 

    \xi_{expval} = \xi_{expval}^I + \xi_{expval}^{II} + \xi_{expval}^{I,II}


where the first two terms involve the contributions from isolated susbsystems (:math:`X\in\{I,II\}`):

.. math::
   :label: fde_magn_expval_11

    \xi_{expval}^{X}  = \int v_{emb}^X \Omega_{ii}^{B,B,X} + \iint w_{emb}^{X,X} \Omega_{ii}^{B,X} \Omega_{jj}^{B,Y}

and the last term involves the coupling between subsystems (:math:`X,Y\in\{I,II\}; X\neq Y`):

.. math::
   :label: fde_magn_expval_12

    \xi_{expval}^{XY} = \iint w_{emb}^{I,II} \Omega_{ii}^{B,X} \Omega_{jj}^{B,Y}

As in the calculations of the property gradient, :math:`v_{emb}^X`, :math:`w_{emb}^{X,X}` and `w_{emb}^{X,Y}` are 
the embedding potential of subsystem X (:math:`X\in\{I,II\}`),
the embedding kernel involving subsystem X (:math:`X\in\{I,II\}`) and the embedding kernel 
describing the coupling between subsystems X and Y (:math:`X,Y\in\{I,II\}; X\neq Y`), respectively.

The first-order perturbed overlap distribution for a subsystem X (:math:`X\in\{I,II\}`) is as in Eq.:ref:`mo_lao_b`, while the second-order derivative
of the overlap distribution is more complex as it involves second-order derivatives of the "direct" and "reorthonormalization" terms
presented in Eq.:ref:`mo_lao_b`, as well as a mixed derivative involving both, "direct" and "reorthonormalization" contributions.





NMR shielding and magnetizability with FDE - sample calculations
================================================================

We will consider a simple system of SeH\ :sub:`2`\ - H\ :sub:`2`\ O, in which SeH\ :sub:`2`\ 
constitutes an active subsystem (subsystem I) and H\ :sub:`2`\ O\ is an environment (subsystem II):

.. image:: h2se_h2o.jpg
 
The computational protocol is different depending on whether the coupling between two
subsystems is included or not (through the term dependent on :math:`w_{emb}^{I,II}`).
To avoid confusion, we will describe these two protocols separately.

We will need two molecular input files. File ``h2se.mol``, corresponding to an active subsystem I, reads:

.. literalinclude:: h2se.mol

and file ``h2o.mol``, standing for a subsystem II ('environment'):

.. literalinclude:: h2o.mol

In this example we will use the four-component relativistic DC Hamiltonian and DFT/LDA with small grid and small basis sets.



NMR shielding calculations without the FDE coupling terms
---------------------------------------------------------
.. _numgrid

Here we demonstrate how to include the first two terms of Eq. :eq:`fde_pg` in the property gradient: 
terms dependent on the embedding potential and the uncoupled embedding kernel on an example of 
the NMR shielding tensor of Se nucleus. 

In this tutorial, we use a coarse grid generated in DIRAC in preceeding calculations
on SeH\ :sub:`2`\ - H\ :sub:`2`\ O 'supermolecule', which was then exported on a ``numerical_grid`` file.
In DIRAC the ``numerical_grid`` file, which can be retrieved from DFT calculations, is divided into 
batches of points. However, FDE module is able to read only the continuous list of grid points,
therefore such ``numerical_grid`` file has to be reformatted.
This can be done by the script which can be found in ``utils/fde-mag`` directory of DIRAC,
``dft_to_fde_grid_convert.py``, which should be called from the directory where the ``numerical_grid`` file is:

.. literalinclude:: grid_convert

This script will overwrite the ``numerical_grid`` and save its original version as ``numerical_grid.ORIG``
(if ``numerical_grid.ORIG`` exists in your working directory, then the script will ask if it should be overwritten).


Step 1: SCF calculation on the subsystem II
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. _step1

The first step is the calculation on subsystem II, the aim of which is
to export its density and electrostatic potential on a text file.

The menu file ``getfrozen.inp`` reads:

.. literalinclude:: getfrozen.inp

and we make sure to get the text file ``GRIDOUT`` (related to the :ref:`GRIDOU` keyword).
We should also remember to supply two files, 
related to keywords :ref:`GRIDOU` (``GRIDOUT``) and :ref:`EXONLY` (``FILEEX``),
which are copies of a numerical grid file:

.. literalinclude:: getfrozen.run

After this step, the ``GRIDOUT`` file will be updated: additional columns on that file
correspond to the electrostatic potential and the density of subsystem II
calculated in each point of a grid.

Step 2: SCF and response calculation on the subsystem I
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. _step2

Now we can use the information about the subsystem II written on the ``GRIDOUT`` file for constructing
the embedding potential, :math:`v_{emb}^{I}`, and embedding kernel, :math:`w_{emb}^{I,I}`.
In order to do that, we use the ``GRIDOUT`` file saved from Step 1 and copy it to a ``FROZEN`` file,
which is given as an argument to the :ref:`FRDENS` keyword.

We use the following input file for the SCF step:

.. literalinclude:: updatefrozen.inp

and the ``pam`` command is:

.. literalinclude:: updatefrozen.run



To calculate the NMR shielding tensor, we use the following input file:

.. literalinclude:: shield_v_w11.inp

in which the keyword :ref:`LAO11` denotes that the contributions from the embedding potential, :math:`v_{emb}^{I}`,
and the uncoupled embedding kernel, :math:`w_{emb}^{I,I}`, will be added to the property gradient.
The corresponding ``pam`` command is:

.. literalinclude:: shield_v_w11.run

This is the advised option. However, as the calculations of the uncoupled embedding kernel terms in the property
gradient may be expensive, one may want to restrict the FDE-LAO contributions to the part dependent on the
embedding potential only. 
In this case, the input file reads:

.. literalinclude:: shield_v.inp

and the corresponding ``pam`` command is:

.. literalinclude:: shield_v.run


Additionally, we may want to save the ``PAMXVC`` and ``TBMO`` files, if we are interested in
visualization of the NMR shielding density and/or magnetically-induced current density.
The calculation on our test molecule give the following results for the isotropic part of the shielding tensor of Se atom:

    ===================================  =========
    terms included in property gradient    value
    ===================================  =========
    None                                 2469.8091
    v                                    2370.3941
    v+w11                                2370.2488
    ===================================  =========

``None`` in the table denotes calculations performed with no "new" FDE-LAO terms in the property gradient (the gauge origin is in the center of mass of SeH2).
Also, in this test case the value obtained with the embedding potential and the uncoupled embedding kernel in the property gradient (``v+w11``) 
is very close to the result of calculations with the embedding potential only (``v``), but for heavier nuclei the difference between these two setups may
be significant.


Calculations with the FDE coupling terms
----------------------------------------

Here we will demonstrate how to include the coupling between two subsystems in the property gradient 
(term dependent on :math:`w_{emb}^{I,II}`).

Step 1: SCF and response calculation on the subsystem II
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the first step we perform the calculations on subsystem II, the aim of which is
to get its density and electrostatic potential written on a text file.
In addition, we will export the perturbed density of subsystem II. 
We will need the DFCOEF file of subsystem II, which can be taken from ``getfrozen.inp`` run with the following input:

.. literalinclude:: getfrozen.inp

with the ``pam`` command:

.. literalinclude:: getfrozen.run

(exactly as done in calculations without the FDE coupling terms presented above). 
From this step we save the ``GRIDOUT.h2o`` file.
In order to get the density perturbed by an external magnetic field in LAO basis, the following input can be used:

.. literalinclude:: getfrozenpert.inp

with the keyword :ref:`EXPPED` (under :ref:`NMR` section) asking the program to save the contributions to the perturbed density on two files ``pertden_direct_lao.FINAL`` and ``pertden_reorth_lao.FINAL``. The corresponding pam command is:

.. literalinclude:: getfrozenpert.run

NOTE: This step has to be run in serial mode for the moment.


Step 2: Response calculation on the subsystem I
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Now we can finally calculate the FDE contributions to the NMR shielding tensor of an active subsystem (H2Se) 
which depend on the embedding kernel including the coupling terms between two subsystems. 
If the SCF step on an active subsystem was not done before, we should run the following input:

.. literalinclude:: updatefrozen.inp

with the ``pam`` command:

.. literalinclude:: updatefrozen.run

The response input is then the following:

.. literalinclude:: shield_v_w11_w12_nocoulomb.inp

with two new keywords: :ref:`LAOFRZ` and :ref:`PERTIM`. The corresponding ``pam`` command is:

.. literalinclude:: shield_v_w11_w12_nocoulomb.run

This input does not include the computationally expensive Coulomb term in the embedding kernel, only the nonadditive
exchange-correlation and kinetic parts of the embedding kernel. 
If the Coulomb term is to be calculated then an additional keyword is required in the input, :ref:`LFCOUL`.
Instead of using the :ref:`LFCOUL` and the :ref:`LAOFRZ` keywords together, it is possible to use a single :ref:`LAO12` keyword,
which represents both contributions to the coupled embedding kernel:

.. literalinclude:: shield_v_w11_w12.inp

with an analogous ``pam`` command:

.. literalinclude:: shield_v_w11_w12.run

This input contains all the necessary contributions to the property gradient.
Additionally, we may want to save the ``PAMXVC`` and ``TBMO`` files, if we are interested in
visualization of the NMR shielding density and/or magnetically-induced current density.
With this setup the calculations on the test system gave the following results:


    ===================================  =========
    terms included in property gradient    value
    ===================================  =========
    v+w11+w12, no Coulomb                2370.2536
    v+w11+w12, ALL terms                 2370.2536
    ===================================  =========

Here the two results are the same within 0.0001 ppm accuracy, as the contributions from the Coulomb embedding contribution to the perturbed Fock matrix 
are very small (can be printed in the output with the print level set at least to 5 under ``**PROPERTIES`` section).



Magnetizability calculations without the FDE coupling terms
-----------------------------------------------------------
In this section we demonstrate how to perform the FDE calculations of the magnetizability tensor, in which the terms dependent on the embedding kernel 
involving the coupling between subsystems are neglected, both in the property gradient and in the expectation value parts.

Analogously to the NMR shielding calculations, the magnetizability calculations are done in few steps. First steps are common for both properties.
In what follows we assume that:
  * the ``numerical_grid`` file is prepared as described in section :ref:`numgrid`,
  * the "Step 1" discussed in section :ref:`step1` is completed, therefore the ``GRIDOUT`` file with the exported charge density of subsystem II ("frozen") is available,
  * the SCF part of the "Step 2" discussed in section :ref:`step2` is completed, therefore the DFCOEF file for subsystem I is available ("DFCOEF.h2se")

Step 2: Response calculation on the subsystem I
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To calculate the magnetizability tensor, we use the following input file:

.. literalinclude:: magn_v_w11.inp

in which the keyword :ref:`LAO11` denotes that the contributions from the embedding potential, :math:`v_{emb}^{I}`,
and the uncoupled embedding kernel, :math:`w_{emb}^{I,I}`, will be added to the property gradient.
The keyword :ref:`EMAFDE` is added to turn on the FDE-LAO terms in the expectation value part of the magnetizbaility tensor,
while the keyword :ref:`MNONKE` excludes the contributions dependent on the non-additive part of the coupled embedding kernel.

The corresponding ``pam`` command is:

.. literalinclude:: magn_v_w11.run

This is the advised option. However, it is also possible to restrict the FDE-LAO contributions (in both, property gradient and in
the expectation value term) to the part dependent on the embedding potential only. 
In this case, the input file reads:

.. literalinclude:: magn_v.inp

and the corresponding ``pam`` command is:

.. literalinclude:: magn_v.run

The calculation on our test molecule give the following results for the isotropic part of the magnetizability tensor of both subsystems:

    ===================================  =========  =========  =========
    terms included in property gradient    H2Se       H2O       H2Se-H2O
    ===================================  =========  =========  =========
    Isolated                             -7.0759     -2.6733   -9.7592
    None                                 -7.2109     -2.6989   −9.9098               
    v                                    -7.1182     -2.6717   −9.7899                   
    v+w11                                -7.1182     -2.6717   -9.7899                      
    supermolecule                                              -10.0294
    ===================================  =========  =========  =========

``None`` in the table denotes calculations performed with no "new" FDE-LAO terms (neither in the property gradient nor in the expectation value part,
the gauge origin is in the center of mass of SeH2-H2O). ``Supermolecule`` value refers to calculations on the full H2Se-H2O system, in other cases the
value shown in the ``H2Se-H2O`` column is a sum of values obtained for the individual subsystems.



Magnetizability calculations with the FDE coupling terms
--------------------------------------------------------







