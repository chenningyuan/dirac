:orphan:

.. _effQED:

===========================================================
Effective QED potentials
===========================================================

QED effects are high-order relativistic effects that cannot be captured in the framework
of the no-pair approximation. A rough explanation of the QED effects that
contribute to the Lamb shift of bound-state electrons is as follows:

-  Vacuum polarization (VP): A charge in space is surrounded by virtual electron-positron pairs
-  Self-energy (SE): The electron drags along its electromagnetic field

Effective QED potential is an approximate way to include QED effects to all-electron relativistic Hamiltonians. In the current implementation,
the following potentials are available: For more information see Ref. :cite:`Sunaga2022`.

Vacuum polarization
================================

1) Uehling potential
------------------------------------------

This the lowest-order term, but the dominant contribution to the vacuum polarization
effect. This is directly derived from the scattering matrix of bound-state electrons :cite:`Uehling1935`.

.. math::
   \varphi_{\text {Ueh. }} (\boldsymbol{x})= & \frac{Z e}{4 \pi \varepsilon_0 r_x} \bar{\lambda} \frac{2 \alpha}{3} \int_0^{\infty} r_y d r_y \rho^{\text {nuc. }}\left(r_y\right) \\ & \times\left[K_0\left(\frac{2}{\bar{\lambda}}\left|r_x-r_y\right|\right)-K_0\left(\frac{2}{\bar{\lambda}}\left|r_x+r_y\right|\right)\right]

where

.. math::

   K_0(x)=\int_1^{\infty} d \zeta e^{-x \zeta}\left(\frac{1}{\zeta^3}+\frac{1}{2 \zeta^5}\right) \sqrt{\zeta^2-1}.

   
Self-energy
================================

Actually, the self-energy effect is represented by a nonlocal matrix potential,
but the following local potentials can be approximately derived or fitted to reproduce
the self-energy correction calculations of few-electron atoms. 

1) Flambaum-Ginges' potential
------------------------------------------------------

The following magnetic (:math:`\varphi_{\text {mag }}`) and high-frequency (:math:`\varphi_{\text {HF}}`) potentials are obtained from
the scattering of a free electron due to the interaction with a classical external potential
(the vertex correction problem). The low-frequency part cannot be taken into account in (:math:`\varphi_{\text {HF}}`)
because of the energy parameter in the original equation. The low-frequency term (:math:`\varphi_{\text {LF}}`) is
artificially obtained, assuming that the low-frequency part should have the range of a 1s orbital of
hydrogen-like atoms :cite:`Flambaum2005`. 

.. math::
   \varphi_{\text {FG }} = \varphi_{\text {mag }} + \varphi_{\text {HF }} + \varphi_{\text {LF}},
   
where

.. math::
   \varphi_{\mathrm{mag}} (\boldsymbol{x})=\frac{\alpha \hbar}{4 \pi m_e c} i {\gamma} \cdot \nabla_x\left[ \Phi(\boldsymbol{x}) \left(K_m\left(\frac{2 r_x}{\bar{\lambda}}\right)-1\right)\right],

.. math::
  \varphi_{\text {HF }}(\boldsymbol{x}, \lambda)=-\frac{\alpha}{\pi} A(Z,\boldsymbol{x})\Phi(\boldsymbol{x}) K_e\left(\frac{2 r_x}{\lambda}\right),

.. math::
   \varphi_{\mathrm{LF}} (\boldsymbol{x})=-\frac{B(Z)}{e} Z^4 \alpha^5 m_e c^2 e^{-Z r_x / a_0}.

The functions :math:`K_m(x)` and :math:`K_e(x)` are as follows:

.. math::
   K_m(x)=\int_1^{\infty} d \zeta \frac{e^{-x \zeta}}{\zeta^2 \sqrt{\zeta^2-1}};

.. math::
   K_e(x)= & \int_1^{\infty} d \zeta \frac{e^{-x \zeta}}{\sqrt{\zeta^2-1}}\left\{-\frac{3}{2}+\frac{1}{\zeta^2}+\left(1-\frac{1}{2 \zeta^2}\right)\right. \\ & \left.\times\left[\ln \left(\zeta^2-1\right)+  4 \ln \left(\frac{1}{Z \alpha}+\frac{1}{2}\right)   \right]\right\} .

:math:`A` and :math:`B` are fitting functions.
:math:`\Phi(\boldsymbol{x})` is the nuclear electrostatic potential,
which is :math:`\frac{Z e}{4 \pi \varepsilon_0 r_x}` in the case of the point-charge model.
The coefficients in :math:`A` and :math:`B` are
adjusted to reproduces the self-energy corrections of the hydrogen-like atoms.
The above equation is correct in the case of the point nucleus, but not the exact form of
the finite-size nuclear charge model.
However, it is confirmed that this approximation works well in the case of the Uehling potential :cite:`Artemyev1997`. 


2) Pyykko-Zhao's potential 
-----------------------------------------------

This is a Gaussian-type model potential where the parameters are fitted to reproduce
the orbital energy and M1 hyperfine splitting within the range  :math:`29 \geq Z \geq 83` :cite:`Pyykko2003`.
      
.. math::
   \varphi_{\mathrm{PZ}}(\boldsymbol{x})=B(Z) e^{-\beta(Z) r_x^2}

:math:`B` and :math:`\beta` are fitting functions.

   
How to use
================================

Effective QED potentials, which are one-electron operators,
should be added to the one electron Hamiltonian, as follows:

.. math::
   \hat{h}_D \rightarrow \hat{h}_D  -e \varphi_{\mathrm{effQED}}

.. math::
   \varphi_{\mathrm{effQED}}  =\sum_A\left(\varphi_A^{\mathrm{SE}}+\varphi_A^{\mathrm{VP}}\right)
   
Example 1: DCG Hamiltonian with Uehling and FG potentials
-------------------------------------------------------------
.. literalinclude:: ueh_fg.inp

Here the reader can find that the matrix part of the magnetic term is "ALPHADOT (:math:`\alpha`)",
although :math:`\gamma` appears in the above expression of the magnetic term. This input is correct
because the :math:`\beta` part in :math:`\gamma=\beta\alpha` is already included
in the code for a technical reason.

Example 2: DCG Hamiltonian with Uehling and PZ potentials
-------------------------------------------------------------
.. literalinclude:: ueh_pz.inp


Note on the application to core-properties
================================

The above two SE potentials give very similar results in the case of valence properties :cite:`Sunaga2022`, but they may give different results in the case of core-properties :cite:`Sunaga2021`.


