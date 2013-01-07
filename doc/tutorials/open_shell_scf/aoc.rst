Average-of-configuration open-shell Hartree-Fock
================================================

This is a short introduction to the theory behind average-of-configuration open-shell Hartree-Fock as implemented in DIRAC. For a more complete description the 
reader may consult chapter 3 of the PhD thesis of Jørn Thyssen [Thyssen2004]_ .

It should first be noted that there is no restricted open-shell Hartree-Fock (ROHF) code in DIRAC. 
The reason is that spin-orbit interaction couples spin and spatial degrees of freedom and make the
formalism much more complicated since one cannot exploit spin symmetry alone for fixing the expansion
coefficients in the reference configuration state function (CSF) which serves as the trial function.

Instead of optimizing the energy for a single open-shell state, we shall optimize the energy for a limited set 
of open-shell states. 


Energy expression
-----------------

Suppose that we have a set of :math:`N_{det}` of Slater determinants, constituting our 
N-particle basis. We next construct and diagonalize a CI matrix in this basis. This gives :math:`N_{det}` 
solutions of the form

.. math::

   \left|\Psi_I\right> = \sum_{P=1}^{N_{det}} \left|\Phi_P\right> C_{PI}

We will now find the set of orbitals which minimizes the average energy

.. math::

   E^{ave}=\frac{1}{N_{det}}\sum_{I=1}^{N_{det}}\left<\Psi_P\left|H\right|\Psi_P\right>

Inserting the expansion of the solutions in terms of Slater determinants and using the fact that the expansion
coefficients :math:`C_{PI}` are  elements of a unitary matrix we obtain

.. math::

  E^{ave}=\frac{1}{N_{det}}\sum_{P=1}^{N_{det}}\sum_{Q=1}^{N_{det}}\left<\Phi_P\left|H\right|\Phi_Q\right>\sum_{I=1}^{N_{det}}C_{PI}^*C_{QI}
         =\frac{1}{N_{det}}\sum_{P=1}^{N_{det}}\left<\Phi_P\left|H\right|\Phi_P\right>

showing that the average can also be taken over the N-particle basis itself.

Introducing open shells and active electrons
--------------------------------------------

The above average energy expression is a functional of the orbitals entering the Slater determinants. We will make a distinction between :

- *inactive orbitals*, present in all Slater determinants, represented by indices :math:`ijkl`
- *active orbitals*, present in some, but not all Slater determinants, represented by indices :math:`uvxy`
- *secondary orbitals*, not present in any Slater determinant, represented by indices :math:`abcd`

We shall also employ indices :math:`pqrs`  for general orbitals.

In order to generate our N-particle basis for averaging we will distribute the orbitals into a number of shells. Each shell :math:`S` is specified by
:math:`M_S` orbital and :math:`N_S` electrons. Inactive and secondary shells have :math:`N_S=M_S` and :math:`N_S=0`, respectively, whereas
active shells have :math:`N_S < M_S`. We generate our N-electron basis by distributing all active electrons in all possible ways within their
respective shells. The total energy can then be written in terms of orbitals rather than Slater determinants as

.. math::

   E^{ave}=\sum_S f_s \left[\sum_{p\in S}h_{pp}+\frac{1}{2}\sum_{T}Q^{T}_{pp}+\frac{1}{2}(a_S-1)Q^S_{pp}\right]

where we have introduced

- the fractional occupation :math:`f_S=\frac{N_S}{M_S}` of shell :math:`S`
- the coupling coefficient :math:`a_S`

  - :math:`a_S=\frac{M_S\left(N_2-1\right)}{N_S\left(M_2-1\right)}` for :math:`f_s\ne 0` 
  - :math:`a_S=1` for :math:`f_s = 0`

- scaled anti-symmetrized two-electron integrals :math:`Q^S_{pq}=f_s\sum_{r\in S}\left[(pq|rr)-(pr|rq)\right]`

Orbital rotations
-----------------

We will use a exponential parametrization for the rotation of orbitals

.. math::

   \tilde{\phi}_p = \sum_q\phi_q\left[exp(-\kappa)\right]_{pq}

where :math:`\kappa` is an anti-Hermitian matrix to ensure unitarity of the transformation.
The exponential parametrization allows for unconstrained optimization (no Lagrange multipliers).
It also allows the easy identification of redundant variational parameters, that is, parameters
whose variation does not change the energy. In this particular case one finds that rotations 
within shells are redundant and the corresponding matrix elements :math:`\kappa_{pq}` can be set 
to zero. 

Gradient elements and off-diagonal blocks of the Fock matrix
------------------------------------------------------------

The generally non-zero elements of the gradient vector are:

- inactive-secondary rotations:
.. math:: \quad g_{ia}=h_{ai}+\sum_T Q_{ai}^T=F^I_{ai}
- active-secondary rotations:
.. math:: \quad g_{ua}=f_U\left[F^I_{au}+(a_U-1)Q^U_{au}\right];\quad u\in U
- inactive-active rotations:
.. math:: \quad g_{iu}=(1-f_U)\left[F^I_{ui}+f_U\alpha_UQ^U_{ui}\right];\quad u\in U
- inter-shell active-active rotations
.. math:: \quad g_{uv}=(f_U-f_V)F^I_{vu}-(a_U-1)f_UQ^U_{vu}-(a_V-1)f_VQ^V_{vu};\quad u\in U,\ v\in V

These gradient elements allow the definition of the off-diagonal elements of the Fock matrix:

- inactive-secondary block:
.. math:: \quad {\cal F}_{ia}=F^I_{ia}
- active-secondary block:
.. math:: \quad {\cal F}_{ua}=F^I_{ua}+(a_U-1)Q^U_{ua};\quad u\in U
- inactive-active block:
.. math:: \quad {\cal F}_{iu}=F^I_{iu}+f_U\alpha_UQ^U_{iu};\quad u\in U
- inter-shell active-active rotations :math:`U\ne V`:
.. math:: \quad {\cal F}_{uv}=F^I_{uv}-\frac{(a_U-1)}{(f_U-f_V)}Q^U_{uv}-\frac{(a_V-1)}{(f_V-f_U)}Q^V_{uv};\quad f_U\ne f_V
.. math:: \quad {\cal F}_{uv}=(a_U-1)Q^U_{uv}-(a_V-1)Q^V_{vu};\quad f_U = f_V

Diagonal blocks of the Fock matrix
----------------------------------

The diagonal blocks of the Fock matrix are *a priori* not related to gradient elements and there is
therefore freedom of choice in their specification. The specific choice will not affect the total
energy, but will affect orbitals energies as well as convergence of the AOC HF calculation.







