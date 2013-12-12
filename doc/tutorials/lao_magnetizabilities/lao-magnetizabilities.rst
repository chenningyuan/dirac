=============================================
Magnetizabilities with London Atomic Orbitals
=============================================

Introduction
------------

In this tutorial we will look at the calculation of (static) magnetizabilities with DIRAC.
For more details about the method, see [Ilias2013]_.

The component :math:`\zeta_{\alpha\beta}` of the static :math:`3\times 3` magnetizability tensor describes
connects component :math:`\alpha` of first-order induced magnetic dipole to component :math:`\beta` of the
external inducing homogeneous magnetic field

.. math::

   m^{(1)}_\alpha = \sum_\beta \zeta_{\alpha\beta} B_{\beta} = \frac{1}{2}\int\left(\mathbf{r}_G\times\mathbf{j}^{B_\beta}\right)_{\alpha}d\tau

The first-order induced magnetic dipole is proportional to the direct product of the position vector :math:`\textbf{r}_G` with 
repect to some arbitrary (gauge) origin :math:`G` and the first-order induced current density :math:`\mathbf{j}^{\mathbf{B}}`.
Although the static induced magnetic dipole is formally independent of the gauge-origin this does not hold true in the finite 
basis approximation, where excruciating slow convergence of the  magnetizability with respect to basis set is typically observed.
London Atomic Orbitals (LAOs), also known as Gauge Including Atomic Orbitals (GIAOs), remove any reference to the arbitrary gauge 
origin by shifting the gauge origin to the centers of individual basis functions and effectively cure both problems.

Example: :math:`NF_3`
---------------------

We shall illustrate these features using the nitrogen trifluoride molecule. 

.. image:: NF3.jpg
   :scale: 50
   :alt: alternate text

The molecular input file `NF3.mol` uses the experimental geometry and automatic symmetry detection.

.. literalinclude:: NF3.mol

Generating the HF wave function
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We first run a Hartree-Fock calculation to generate orbitals and corresponding energies
using `scf.inp`

.. literalinclude:: scf.inp

and the command::

   pam --inp=scf --mol=NF3 --outcmo

You may notice in the output that DIRAC will center and rotate the molecule.
It detects the full :math:`C_{3v}` symmetry, but will run the calculation in 
the lower :math:`C_s` symmetry, with reflection in the *xy* -plane.

Magnetizabilities: first attempt
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We first calculate the magnetizability using the input file `cgo.inp`

.. literalinclude:: cgo.inp

where the gauge origin has been set to the center of mass using the :ref:`.USECM` keyword.
Notice that we use :ref:`.RKBIMP` to convert our molecular coefficients from restricted
to unrestricted kinetic balance (RKB :math:`\rightarrow` UKB), the former employed for the
generation of orbitals, and the latter employed for the response calculations, indicated by
:ref:`.URKBAL`. This corresponds to the use of `simple magnetic balance <../simple_magnetic_balance/tutorial.html>`_.
The calculation is run using::

   pam --inp=cgo --mol=NF3 --incmo

and gives the total magnetizability tensor

.. literalinclude:: cgo_NF3_total

here reported in atomic units :math:`e^2a_0^2/m_e`, corresponding to 7.89104 :math:`\cdot 10^{-29}` J/T.

Magnetizabilities using LAOs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We now activate LAOs using the input file `lao.inp`

.. literalinclude:: lao.inp





