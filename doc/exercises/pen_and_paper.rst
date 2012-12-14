.. _day_1_exercises:

The course comprises two pen-and-paper exercises which are based on materials of the morning-session lectures of Day 1.

Dirac's relation
================

A relation that is often exploited throughout the course, e.g. to achieve a separation of spin-dependent and
independent terms in the Dirac equation (see lecture on Day 2), is Dirac's relation, which can be
written for two arbitrary vector operators :math:`\vec{u}` and :math:`\vec{v}` as:

.. math::

  (\vec{\sigma} \cdot \vec{u})(\vec{\sigma} \cdot \vec{v}) = \vec{u} \cdot \vec{v} I_{2} + i \vec{\sigma} \cdot (\vec{u} \times \vec{v})

where :math:`\vec{\sigma}` are the Pauli spin matrices and :math:`I_{2}` is a :math:`2 \times 2` unit matrix. Note that :math:`\vec{u}` and :math:`\vec{v}` do not necessarily commute.

* **Problem 1**: Verify the relation expanding the left-hand side of Dirac's relation. Hint: use the relations for the Pauli spin matrices:

.. math::

  \vec{\sigma} = (\sigma_x, \sigma_y, \sigma_z)

  \sigma_x^2 = \sigma_y^2 = \sigma_z^2 = I_{2}

  \sigma_i \sigma_j = \delta_{ij} I_{2} + i \sum_{k=1}^{3} \epsilon_{ijk} \sigma_k


* **Problem 2**: derive a final expression inserting for :math:`\vec{u} = \vec{v}` the kinematical momentum operator :math:`\vec{\pi} = \vec{p} - \frac{q_e}{c}\vec{A}` (where :math:`\vec{A}`  is an external electromagnetic vector potential). Hint:

.. math::

  [\pi_i,\pi_j] = \frac{q}{c} i \hbar \sum\limits_{k=1}^{3} \epsilon_{ijk} B_k


Two-component Pauli equation (0th order Pauli equation)
=======================================================

From the one-electron Dirac equation in external magnetic fields we may derive the Pauli equation (also known as 0th order Pauli
equation) by considering the non-relativistic limit :math:`c \rightarrow \infty`.

.. math::

  \left[ \frac{\vec{p}^2}{2m_e} + \frac{q_e^2 \vec{A}^2}{2m_e c^2} - \frac{q_e}{2m_ec}(\bf{l} + 2 \bf{s})
  \cdot \vec{B} + V \right] \Psi^L = i \hbar  \frac{\partial}{\partial t} \Psi^L

* **Problem 3**: Derive the above equation starting from the one-electron Dirac equation in external magnetic fields written in two-spinor form:

.. math::

  i \hbar \frac{\partial}{\partial t} \left( \begin{array}{c}
  \Psi^L \\
  \Psi^S \end{array} \right) = c  \left( \begin{array}{c}
  (\vec{\sigma} \cdot \vec{\pi}) \Psi^S \\
  (\vec{\sigma} \cdot \vec{\pi}) \Psi^L \end{array} \right)
  + m_ec^2  \left( \begin{array}{c}
   \Psi^L \\
  -\Psi^S \end{array} \right) + V  \left( \begin{array}{c}
  \Psi^L \\
  \Psi^S \end{array} \right)

The following hints may be useful:

1. shift the energy limit to the non-relativistic limit (E = 0 rather than E = :math:`m_ec^2`):

.. math::

  V \rightarrow V - m_ec^2

2. elimininate the small-component :math:`\Psi^S`  from the upper component using the kinetic balance condition:

.. math::

  \Psi^S \approx \frac{\vec{\sigma} \cdot \vec{\pi}}{2m_ec^2} \Psi^L


3. use :math:`\vec{A} \cdot \vec{p} = \vec{p} \cdot \vec{A} - (\vec{p} \cdot \vec{A})`

4. remember for the Coulomb gauge:

.. math::

 div(\vec{A}) = 0

5. for a constant and homogeneous magnetic field :math:`\bf{B} = curl(\bf{A})`, we may write :math:`\bf{A} = \frac{1}{2}(\bf{B} \times \bf{r})`

* **Problem 4**: define the *Bohr magneton* and the *gyromagnetic ratio g* of the electron (Lande factor) according to the Pauli Hamiltonian.


Literature and further reading
==============================

* K. Dyall and K. Fœgri, Introduction to Relativistic Qunatum Chemistry, Chapter 4.

* M. Reiher and M. Wolf, Relativistic Quantum Chemistry, Chapter 5.

