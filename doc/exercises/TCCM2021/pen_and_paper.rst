:orphan:

.. _TCCM_pp_exercises:

In the text full units are used for clarity, but in practice one can employ SI-based atomic units, setting :math:`\hbar=m_{e}=e=4\pi\varepsilon_{0}=1`.

Pauli spin matrices
===================

The Pauli spin matrices 

.. math::

   \sigma_{x}=\left[\begin{array}{cc}
   0 & 1\\
   1 & 0
   \end{array}\right];\quad\sigma_{y}=\left[\begin{array}{cc}
   0 & -i\\
   i & 0
   \end{array}\right];\quad\sigma_{z}=\left[\begin{array}{cc}
   1 & 0\\
   0 & -1
   \end{array}\right]

are representation matrices of of the operator :math:`2\mathbf{s}` in the basis :math:`\left\{ \left|\alpha\right\rangle ,\left|\beta\right\rangle \right\}` of spin-1/2 functions (``spin-up'' and ``spin-down'', respectively), that is :math:`2\mathbf{s}\rightarrow\hbar\boldsymbol{\sigma}`. These spin functions obey the following relations

.. math::

   \begin{array}{lclclcl}
   \hat{s}_{z}\left|\alpha\right\rangle  & = & +\frac{1}{2}\hbar\left|\alpha\right\rangle  &  &    \hat{s}_{z}\left|\beta\right\rangle  & = & -\frac{1}{2}\hbar\left|\beta\right\rangle \\
   \hat{s}_{+}\left|\alpha\right\rangle  & = & 0 &  & \hat{s}_{+}\left|\beta\right\rangle  & = &  \hbar\left|\alpha\right\rangle \\
   \hat{s}_{-}\left|\alpha\right\rangle  & = & \hbar\left|\beta\right\rangle  &  &   \hat{s}_{-}\left|\beta\right\rangle  & = & 0
   \end{array}

where appears the ladder operators :math:`s_{\pm}=s_{x}\pm is_{y}`. The spin-1/2 functions are furthermore orthonormal. Let us see how the Pauli :math:`\sigma_{z}` matrix is generated:

 * Element (1,1): :math:`\left\langle \alpha\left|2s_{z}\right|\alpha\right\rangle /\hbar=\left\langle \alpha\left|\right.\alpha\right\rangle =1`
 * Element (1,2): :math:`\left\langle \alpha\left|2s_{z}\right|\beta\right\rangle /\hbar=\left\langle \alpha\left|\right.\beta\right\rangle =0`
 * Element (2,1): :math:`\left\langle \beta\left|2s_{z}\right|\alpha\right\rangle \hbar=\left\langle \beta\left|\right.\alpha\right\rangle =0`
 * Element (2,2): :math:`\left\langle \beta\left|2s_{z}\right|\beta\right\rangle /\hbar=-\left\langle \beta\left|\right.\beta\right\rangle =-1`

        
 1. In the same manner construct representation matrices for the ladder operators :math:`s_{+}` and :math:`s_{-}` and from this :math:`\sigma_{x}` and :math:`\sigma_{y}`.
 2. Demonstrate, for general vector operators :math:`\boldsymbol{A}` and :math:`\boldsymbol{B}`,  the Dirac identity

 .. math ::

         \left(\boldsymbol{\sigma}\cdot\boldsymbol{A}\right)\left(\boldsymbol{\sigma}
         \cdot\boldsymbol{B}\right)=\boldsymbol{A}\cdot\boldsymbol{B}+i\boldsymbol{\sigma}
         \cdot\left(\boldsymbol{A}\times\boldsymbol{B}\right)

 3. Use the Dirac identity to calculate

   * :math:`\left(\boldsymbol{\sigma}\cdot\hat{\boldsymbol{p}}\right)\left(\boldsymbol{\sigma}\cdot\hat{\boldsymbol{p}}\right)`
   * :math:`\left(\boldsymbol{\sigma}\cdot\hat{\boldsymbol{\pi}}\right)\left(\boldsymbol{\sigma}\cdot\hat{\boldsymbol{\pi}}\right)`, where :math:`\hat{\boldsymbol{\pi}}` is mechanical momentum :math:`\hat{\boldsymbol{\pi}}=\hat{\boldsymbol{p}}+e\boldsymbol{A}` and :math:`\boldsymbol{A}` is a vector potential. Show how this expression is simplified in Coulomb gauge: :math:`\boldsymbol{\nabla}\cdot\boldsymbol{A}=0`.

   * Show that :math:`\left(\boldsymbol{\sigma}\cdot\hat{\boldsymbol{p}}\right)V_{eN}\left(\boldsymbol{\sigma}\cdot\hat{\boldsymbol{p}}\right)=\boldsymbol{p}V_{eN}\cdot\boldsymbol{p}+\hbar\left(\frac{Ze^{2}}{2\pi\varepsilon_{0}r^{3}}\right)\hat{\boldsymbol{s}}\cdot\hat{\boldsymbol{\ell}}`, where :math:`V_{eN}` is the electron-nucleus interaction :math:`V_{en}=-\frac{Ze^{2}}{4\pi\varepsilon_{0}r}` in an atom


 
   


