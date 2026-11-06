:orphan:
 

================================
MP2 Natural Orbital tutorial
================================

Introduction
============

The coupled-cluster (CC) method is quite expansive. In the case of the CCSD method, it is of :math:`O(n_o^2n_v^4)` with :math:`n_o` is the number of occupied orbitals and :math:`n_v` the number of virtual orbitals. A way to reduce the scaling of the CC methods, by a prefactor, is Natural Orbitals (NO). The main point is to reduce the number of virtual orbitals by compressing the virtual space into a smaller more meaningful one, by diagonalizing the virtual-virtual block of the density matrix. The diagonalizing procedure generates eigenpairs of new virtual orbitals (eigenvectors) and occupation numbers (eigenvectors). Then, these new virtual orbitals are truncated by only keeping orbitals whose occupation number is higher than a threshold. In the MP2NO framework, the density matrix to diagonalize is the MP2 density matrix :cite:`Davidson1972`.


MP2NO on HF system
==================

MP2NO has been used in Dirac here :cite:`Yuan2022`. Running a NO calculation is done in two steps:

1. Generate the NO that will be used in step 2. The orbitals will be saved in the `NOs_MO` file.
2. Put the `NOs_MO` file and read the NO from step 1 and compute the desired properties with them

With the following `hf` system:

.. literalinclude:: hf.xyz

We can do step 1 with this input file:

.. literalinclude:: hf_nr.inp

The command to launch step 1 is, please note the `--outcmo` option:

.. code-block:: bash

   $DIRACPAM <DIRAC_OPTIONS> --inp hf_nr.inp --mol hf.xyz --outcmo

The truncation is done by chosing the threshold value under the `.MP2NO` keyword. Afterwards, you need to check the number of virtual orbitals inside the `NOs_MO` file from the output of step 1 that should look like:

.. code-block:: bash

     MP2 occupation number (o.n.) threshold (t) :   0.1000D-02
     Number of MP2 FVNOs retained (o.n. > t)    :            9

     Sum of MP2 FVNOs o.n., retained NOs (a)    :   0.3874D-01
     Sum of MP2 FVNOs o.n., all NOs (b)         :   0.3980D-01
     Percentage of occupation retained (100*a/b):   0.9734D+02


It will be useful for the input file for step 2, for computing the CCSD energy, dipole moment and electric-field-gradient:

.. literalinclude:: property_nr.inp

The command to launch step 2 is (please note the `--incmo` option):

.. code-block:: bash

   $DIRACPAM <DIRAC_OPTIONS> --inp property_nr.inp --mol hf.xyz --incmo
