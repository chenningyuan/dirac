:orphan:


.. _one_electron_operators:

One-electron operators
======================

Syntax for the specification of one-electron operators
------------------------------------------------------

A general one-electron operator in 4-component calculations is generated from linear combinations of operators having the basic form:

.. math::

  \hat{P} = f \, M_{4 \times 4} \, \hat{\Omega}

where :math:`f` is a scalar factor, :math:`\hat{\Omega}` is a scalar operator, and :math:`M_{4 \times 4}` is one of the following :math:`4 \times 4` matrices:

.. math::

  & I_{4 \times 4}, \, \beta, \, \gamma_5, \, \beta\gamma_5,

  & i\alpha_x, \, i\alpha_y, \, i\alpha_z, \, i\beta\alpha_x, \, i\beta\alpha_y, \, i\beta\alpha_z,

  & i\Sigma_x, \, i\Sigma_y, \, i\Sigma_z, \, i\beta\Sigma_x, \, i\beta\Sigma_y, \, i\beta\Sigma_z.

One thing to notice is that the Dirac matrices :math:`I_{4 \times 4}`, :math:`\beta`, :math:`\gamma_5`, and :math:`\beta\gamma_5` are time-reversal symmetric operators, and they are implemented with such symmetry in DIRAC. On the other hand, the operators :math:`\boldsymbol{\alpha}`, :math:`\boldsymbol{\Sigma}`, :math:`\boldsymbol{\beta\alpha}`, and :math:`\boldsymbol{\beta\Sigma}` are time-reversal antisymmetric. This information is saved by the program. In addition, these time-reversal antisymmetric Dirac matrices (:math:`\boldsymbol{\alpha}`, :math:`\boldsymbol{\Sigma}`, :math:`\boldsymbol{\beta\alpha}`, and :math:`\boldsymbol{\beta\Sigma}`) are multiplied by the imaginary unit :math:`i=\sqrt{-1}` in order to make them time-reversal symmetric and hence fit into the quaternion symmetry scheme of DIRAC (see :cite:`Saue1999` and :cite:`Salek2005` for more information).

More details about the implemented operators are provided in Section `Quaternion algebra and time-reversal symmetry`_, where the connection between the Dirac matrices and the quaternion algebra is given.


Operator specification
^^^^^^^^^^^^^^^^^^^^^^

Operators are specified by the keyword :ref:`HAMILTONIAN_.OPERATOR` with the following arguments::

  .OPERATOR
   'operator name'
   operator type keyword
   scalar operator label for each component
   FACTORS
   factors for each component
   CMULT
   COMFACTOR
   common factor for all components

Note that the arguments following the keyword :ref:`HAMILTONIAN_.OPERATOR` must start with a blank. The arguments are optional, except for the operator label.
Component factors as well as the common factor are all one if not specified.


.. _Operator type keywords:

Operator type keywords
^^^^^^^^^^^^^^^^^^^^^^

There are 22 basic operator types implemented in DIRAC. They are listed in the following Table:

===========  ==================================================================================================  ===============  ======================================
**Keyword**  **Operator form**                                                                                   **Nr. factors**  **Time-reversal symmetry (without i)**
===========  ==================================================================================================  ===============  ======================================
DIAGONAL     :math:`f I_{4 \times 4} \hat{\Omega}`                                                               1                Symmetric
XALPHA       :math:`f i\alpha_x \hat{\Omega}`                                                                    1                Antisymmetric
YALPHA       :math:`f i\alpha_y \hat{\Omega}`                                                                    1                Antisymmetric
ZALPHA       :math:`f i\alpha_z \hat{\Omega}`                                                                    1                Antisymmetric
XAVECTOR     :math:`f_1 i\alpha_y \hat{\Omega}_z - f_2 i\alpha_z \hat{\Omega}_y`                                 2                Antisymmetric
YAVECTOR     :math:`f_1 i\alpha_z \hat{\Omega}_x - f_2 i\alpha_x \hat{\Omega}_z`                                 2                Antisymmetric
ZAVECTOR     :math:`f_1 i\alpha_x \hat{\Omega}_y - f_2 i\alpha_y \hat{\Omega}_x`                                 2                Antisymmetric
ALPHADOT     :math:`f_1 i\alpha_x \hat{\Omega}_x + f_2 i\alpha_y \hat{\Omega}_y + f_3 i\alpha_z \hat{\Omega}_z`  3                Antisymmetric
GAMMA5       :math:`f \gamma_5 \hat{\Omega}`                                                                     1                Symmetric
XSIGMA       :math:`f i\Sigma_x \hat{\Omega}`                                                                    1                Antisymmetric
YSIGMA       :math:`f i\Sigma_y \hat{\Omega}`                                                                    1                Antisymmetric
ZSIGMA       :math:`f i\Sigma_z \hat{\Omega}`                                                                    1                Antisymmetric
XBETASIG     :math:`f i\beta \Sigma_x \hat{\Omega}`                                                              1                Antisymmetric
YBETASIG     :math:`f i\beta \Sigma_y \hat{\Omega}`                                                              1                Antisymmetric
ZBETASIG     :math:`f i\beta \Sigma_z \hat{\Omega}`                                                              1                Antisymmetric
XBETAALP     :math:`f i\beta \alpha_x \hat{\Omega}`                                                              1                Antisymmetric
YBETAALP     :math:`f i\beta \alpha_y \hat{\Omega}`                                                              1                Antisymmetric
ZBETAALP     :math:`f i\beta \alpha_z \hat{\Omega}`                                                              1                Antisymmetric
BETA         :math:`f \beta \hat{\Omega}`                                                                        1                Symmetric
SIGMADOT     :math:`f_1 i\Sigma_x \hat{\Omega}_x + f_2 i\Sigma_y \hat{\Omega}_y + f_3 i\Sigma_z \hat{\Omega}_z`  1                Antisymmetric
BETAGAMMA5   :math:`f \beta \gamma_5 \hat{\Omega}`                                                               1                Symmetric
iBETAGAMMA5  :math:`f \beta \gamma_5 \hat{\Omega}`                                                               1                Antisymmetric
===========  ==================================================================================================  ===============  ======================================

The last column in the previous table reports the time-reversal symmetry of the operators without imaginary phase.

Note for advanced users: iBETAGAMMA5 corresponds to the time-reversal symmetric operator :math:`f \beta \gamma_5 \hat{\Omega}`, but has been assigned a time-reversal antisymmetry. This operator can be used, for example, for calculations of the molecular enhancement factors of the electron electric dipole moments. For more details, please refer to Section `Quaternion algebra and time-reversal symmetry`_.


.. _One-electron scalar operator labels:

One-electron scalar operator labels
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

+-------------+----------------------------------+----------------+----------------------------------------------------------------------------------------+
| **Operator**| **Integral description**         | **Components** | **Operators**                                                                          |
| **label**   |                                  |                |                                                                                        |
+=============+==================================+================+========================================================================================+
| KINENER     | NR kinetic energy                | KINENER        | :math:`\hat{\Omega}_1 = -\frac{1}{2}\hat{\nabla}^2`                                    |
+-------------+----------------------------------+----------------+----------------------------------------------------------------------------------------+
| MOLFIELD    | Nuclear attraction               | MOLFIELD       | :math:`\hat{\Omega}_1 = \sum_K Z_K |\boldsymbol{r}-\boldsymbol{R}_K|^{-1}`             |
+-------------+----------------------------------+----------------+----------------------------------------------------------------------------------------+
| OVERLAP     | Overlap                          | OVERLAP        | :math:`\hat{\Omega}_1 = 1`                                                             |
+-------------+----------------------------------+----------------+----------------------------------------------------------------------------------------+
| FC          | One-electron Fermi contact       | FC NAMab       | :math:`\hat{\Omega}_1 = \frac{4 \pi g_{e}}{3} \delta(\boldsymbol{r}-\boldsymbol{R}_K)` |
+-------------+----------------------------------+----------------+----------------------------------------------------------------------------------------+
| ANGMOM      | Orbital angular momentum         | XANGMOM        | :math:`\hat{\Omega}_1 = -[(\boldsymbol{r}-\boldsymbol{R}_{CM})\times\hat{\nabla}]_x`   |
|             | around CM                        +----------------+----------------------------------------------------------------------------------------+
|             |                                  | YANGMOM        | :math:`\hat{\Omega}_2 = -[(\boldsymbol{r}-\boldsymbol{R}_{CM})\times\hat{\nabla}]_y`   |
|             |                                  +----------------+----------------------------------------------------------------------------------------+
|             |                                  | ZANGMOM        | :math:`\hat{\Omega}_3 = -[(\boldsymbol{r}-\boldsymbol{R}_{CM})\times\hat{\nabla}]_z`   |
+-------------+----------------------------------+----------------+----------------------------------------------------------------------------------------+
| DIPLEN      | Dipole length                    | XDIPLEN        | :math:`\hat{\Omega}_1 = x`                                                             |
|             |                                  +----------------+----------------------------------------------------------------------------------------+
|             |                                  | YDIPLEN        | :math:`\hat{\Omega}_2 = y`                                                             |
|             |                                  +----------------+----------------------------------------------------------------------------------------+
|             |                                  | ZDIPLEN        | :math:`\hat{\Omega}_3 = z`                                                             |
+-------------+----------------------------------+----------------+----------------------------------------------------------------------------------------+
| DIPVEL      | Dipole velocity                  | XDIPVEL        | :math:`\hat{\Omega}_1 = -\frac{\partial}{\partial x}`                                  |
|             |                                  +----------------+----------------------------------------------------------------------------------------+
|             |                                  | YDIPVEL        | :math:`\hat{\Omega}_2 = -\frac{\partial}{\partial y}`                                  |
|             |                                  +----------------+----------------------------------------------------------------------------------------+
|             |                                  | ZDIPVEL        | :math:`\hat{\Omega}_3 = -\frac{\partial}{\partial z}`                                  |
+-------------+----------------------------------+----------------+----------------------------------------------------------------------------------------+
| SPNORB      | Spatial spin-orbit               | X1SPNORB       |                                                                                        |
|             |                                  +----------------+----------------------------------------------------------------------------------------+
|             |                                  | Y1SPNORB       |                                                                                        |
|             |                                  +----------------+----------------------------------------------------------------------------------------+
|             |                                  | Z1SPNORB       |                                                                                        |
+-------------+----------------------------------+----------------+----------------------------------------------------------------------------------------+
| QUADRUP     | Quadrupole moment                | XXQUADRU       | :math:`\hat{\Omega}_1 = \frac{1}{4} (x^2-r^2)`                                         |
|             |                                  +----------------+----------------------------------------------------------------------------------------+
|             |                                  | XYQUADRU       | :math:`\hat{\Omega}_2 = \frac{1}{4} xy`                                                |
|             |                                  +----------------+----------------------------------------------------------------------------------------+
|             |                                  | XZQUADRU       | :math:`\hat{\Omega}_3 = \frac{1}{4} xz`                                                |
|             |                                  +----------------+----------------------------------------------------------------------------------------+
|             |                                  | YYQUADRU       | :math:`\hat{\Omega}_4 = \frac{1}{4} (y^2-r^2)`                                         |
|             |                                  +----------------+----------------------------------------------------------------------------------------+
|             |                                  | YZQUADRU       | :math:`\hat{\Omega}_5 = \frac{1}{4} yz`                                                |
|             |                                  +----------------+----------------------------------------------------------------------------------------+
|             |                                  | ZZQUADRU       | :math:`\hat{\Omega}_6 = \frac{1}{4} (z^2-r^2)`                                         |
+-------------+----------------------------------+----------------+----------------------------------------------------------------------------------------+
| SECMOM      | Second moment                    | XXSECMOM       | :math:`\hat{\Omega}_1 = - x^2`                                                         |
|             |                                  +----------------+----------------------------------------------------------------------------------------+
|             |                                  | XYSECMOM       | :math:`\hat{\Omega}_2 = - xy`                                                          |
|             |                                  +----------------+----------------------------------------------------------------------------------------+
|             |                                  | XZSECMOM       | :math:`\hat{\Omega}_3 = - xz`                                                          |
|             |                                  +----------------+----------------------------------------------------------------------------------------+
|             |                                  | YYSECMOM       | :math:`\hat{\Omega}_4 = - y^2`                                                         |
|             |                                  +----------------+----------------------------------------------------------------------------------------+
|             |                                  | YZSECMOM       | :math:`\hat{\Omega}_5 = - yz`                                                          |
|             |                                  +----------------+----------------------------------------------------------------------------------------+
|             |                                  | ZZSECMOM       | :math:`\hat{\Omega}_6 = - z^2`                                                         |
+-------------+----------------------------------+----------------+----------------------------------------------------------------------------------------+
| THETA       | Traceless quadrupole             | XXTHETA        | :math:`\hat{\Omega}_1 = - \frac{3}{2} (x^2-\frac{1}{3}r^2)`                            |
|             |                                  +----------------+----------------------------------------------------------------------------------------+
|             |                                  | XYTHETA        | :math:`\hat{\Omega}_2 = - \frac{3}{2} xy`                                              |
|             |                                  +----------------+----------------------------------------------------------------------------------------+
|             |                                  | XZTHETA        | :math:`\hat{\Omega}_3 = - \frac{3}{2} xz`                                              |
|             |                                  +----------------+----------------------------------------------------------------------------------------+
|             |                                  | YYTHETA        | :math:`\hat{\Omega}_4 = - \frac{3}{2} (y^2-\frac{1}{3}r^2)`                            |
|             |                                  +----------------+----------------------------------------------------------------------------------------+
|             |                                  | YZTHETA        | :math:`\hat{\Omega}_5 = - \frac{3}{2} yz`                                              |
|             |                                  +----------------+----------------------------------------------------------------------------------------+
|             |                                  | ZZTHETA        | :math:`\hat{\Omega}_6 = - \frac{3}{2} (z^2-\frac{1}{3}r^2)`                            |
+-------------+----------------------------------+----------------+----------------------------------------------------------------------------------------+

Note that in the one-electron Fermi-contact integrals :math:`g_e` is the electronic g-factor (:math:`g_e = 2` unless QED effects are taken into account), 'NAM' are the three first letters in the name of atom :math:`K` as given in the .mol file, and 'ab' stands for the number of the symmetry-adapted nucleus.


Other one-electron operator labels are:

=========== ====================================================================================================================
**Keyword** **Description**
=========== ====================================================================================================================
ANGLON      Angular momentum around the nuclei
CARMOM      Cartesian moments integrals, symmetric integrals, (l + 1)(l + 2)/2 components ( l = i + j + k)
CM1         First order magnetic field derivatives of electric field
CM2         Second order magnetic field derivatives of electric field
DARWIN      Darwin type integrals
DIASUS      Angular London orbital contribution to diamagnetic susceptibility
DSO         Diamagnetic spin-orbit integrals
DSUSCGO     Diamagnetic susceptibility with common gauge origin
DSUSLH      Angular London orbital contribution to diamagnetic susceptibility
DSUSNOL     Diamagnetic susceptibility without London contribution
ELFGRDC     Electric field gradient at the individual nuclei, cartesian
ELFGRDS     Electric field gradient at the individual nuclei, spherical
EXPIKR      Cosine and sine integrals
HBDO        Half B-differentiated overlap matrix
HDO         Half-derivative overlap integrals
HDOBR       Ket-differentiation of HDO-integrals with respect to magnetic field
LONMOM      London orbital contribution to angular momentum
MAGMOM      One-electron contributions to magnetic moment
MASSVEL     Mass velocity integrals
NEFIELD     Electric field at the individual nuclei
NSTCGO      Nuclear shielding integrals with common gauge origin
NUCPOT      Potential energy of the interaction of electrons with individual nuclei, divided by the nuclear charge
NUCSHI      Nuclear shielding tensor integrals
NUCSLO      London orbital contribution to nuclear shielding tensor integrals
NUCSNLO     Nuclear shielding integrals without London orbital contribution
PSO         Paramagnetic spin-orbit integrals
S1MAG       Second order contribution from overlap matrix to magnetic properties
S1MAGL      Bra-differentiation of overlap matrix with respect to magnetic field
S1MAGR      Ket-differentiation of overlap matrix with respect to magnetic field
SDFC        Spin-dipole + Fermi contact integrals
SOLVENT     Electronic solvent integrals
SPHMOM      Spherical moments integrals (real combinations), symmetric integrals, (2l + 1) components ( m = +0, -1, +1, ..., +l)
SPIN-DI     Spin-dipole integrals
SQHDO       Half-derivative overlap integrals not to be antisymmetrized
SQHDOR      Half-derivative overlap integrals not to be anti-symmetrized
SQOVLAP     Second order derivatives overlap integrals
=========== ====================================================================================================================


Examples of using various operators
-----------------------------------

We give here several concrete examples on how to construct operators for various properties.


Kinetic part of the Dirac Hamiltonian
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The kinetic part of the electronic Dirac Hamiltonian is given by: :math:`\hat{H}^K=c\boldsymbol{\alpha}\cdot\boldsymbol{p}`. This operator may be specified by::

  .OPERATOR
   'Kin energy'
   ALPHADOT
   XDIPVEL
   YDIPVEL
   ZDIPVEL
   COMFACTOR
   137.03599926085

where 137.03599926085 is the CODATA22 value for :math:`c`.

The speed of light :math:`c` is an important parameter in relativistic theory, but its explicit value in atomic units is not necessarily remembered. A simpler way to specify the kinetic energy operator is therefore::

  .OPERATOR
   'Kin energy'
   ALPHADOT
   XDIPVEL
   YDIPVEL
   ZDIPVEL
   CMULT
   COMFACTOR
   1

where the keyword ``CMULT`` assures multiplication of the common factor 1 by :math:`c`. This option has the further advantage that ``CMULT`` follows any user-specified modification of the speed of light, as provided by :ref:`GENERAL_.CVALUE`.


Magnetic interactions
^^^^^^^^^^^^^^^^^^^^^

Another example is the operator that describes the electromagnetic interaction between electrons and a uniform external magnetic field: :math:`\hat{H}^B=\boldsymbol{B}\cdot\left(-\frac{c}{2}\boldsymbol{\alpha}\times\boldsymbol{r}\right)`.

The operator to use in this case is :math:`-\frac{c}{2}\boldsymbol{\alpha}\times\boldsymbol{r}`, and in DIRAC each component is calculated individually; for example, the :math:`x`-component of the operator is :math:`-\frac{c}{2}\left(\alpha_y z - \alpha_z y\right)`. This operator can be requested using::

  .OPERATOR
   'B_x'
   XAVECTOR
   ZDIPLEN
   YDIPLEN
   CMULT
   COMFACTOR
   -0.5

The program will assume all operators to be Hermitian. The operator ``XAVECTOR`` is defined as :math:`f_1 i\alpha_y \hat{\Omega}_z-f_2 i\alpha_z \hat{\Omega}_y` (check `Operator type keywords`_), and in this particular case the scalar operators are :math:`\hat{\Omega}_z=z` and :math:`\hat{\Omega}_y=y` (see ``ZDIPLEN`` and ``YDIPLEN`` in `One-electron scalar operator labels`_). Additionally, :math:`f_1=f_2=-\frac{c}{2}`. Then, DIRAC will assume the operator to be :math:`-\frac{c}{2}\left(i\alpha_y z - i\alpha_z y\right)`.

The imaginary phase :math:`i` makes the operator time-reversal symmetric.


Diagonal operators
^^^^^^^^^^^^^^^^^^

If no other argument is given, the program assumes the operator to be diagonal and expects the operator name to be the component label, for instance::

  .OPERATOR
   OVERLAP


Dipole moment as finite field perturbation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Another example is the finite perturbation calculation with the :math:`\hat{z}` dipole length operator added to the Hamiltonian (don't forget to decrease the symmetry of your system):

.. math::

    \hat{H} = \hat{H}_0 + 0.01 \cdot \hat{z}

::

  .OPERATOR
   ZDIPLEN
   COMFACTOR
   0.01


Fermi-contact integrals
^^^^^^^^^^^^^^^^^^^^^^^

Here is an example where the Fermi-contact (FC) integrals for a certain nucleus :math:`K`
are added to the Hamiltonian in a finite-field calculation.

Let's assume you are looking at a PbX dimer (order in the .mol file: 1. Pb, 2. X) and you want to add the FC integrals for the Pb nucleus to the Dirac-Coulomb :ref:`**HAMILTONIAN` as a perturbation with a given field strength of :math:`10^{-9}`.

The FC integrals are called using the label 'FC NAMab' (see `One-electron scalar operator labels`_), and these operators can be requested using::

  **HAMILTONIAN
  .OPERATOR
   'Density at nucleus'
   DIAGONAL
   'FC Pb 01'
   FACTORS
   0.000000001
   
**Important note:** The values obtained after the fit of the finite-field energies need to be scaled by :math:`\frac{3}{4 \pi g_{e}} = \frac{1}{8.3872954891254192}` in order to get the raw density :math:`\delta(\boldsymbol{r}-\boldsymbol{R}_\text{Pb})`.

The next example shows how to calculate the electronic density at the nucleus as the expectation value :math:`\langle 0 \vert \delta(\boldsymbol{r}-\boldsymbol{R}_K) \vert 0 \rangle` for a Dirac-Coulomb Hartree-Fock wave function including a decomposition of the molecular orbital contributions to the density::

  **DIRAC
  .WAVE FUNCTION
  .PROPERTIES
  **HAMILTONIAN
  **WAVE FUNCTION
  .SCF
  **PROPERTIES
  .RHONUC
  *EXPECTATION
  .ORBANA
  *END OF


Cartesian moment expectation value
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In the following example I am calculating a cartesian moment expectation value :math:`\langle 0 \vert x^1 y^2 z^3 \vert 0 \rangle` for a Levy-Leblond HF wave function::

  **DIRAC
  .WAVE FUNCTION
  .PROPERTIES
  **HAMILTONIAN
  .LEVY-LEBLOND
  **WAVE FUNCTION
  .SCF
  **PROPERTIES
  *EXPECTATION
  .OPERATOR
   CM010203
  *END OF


.. _Quaternion algebra and time-reversal symmetry:

Quaternion algebra and time-reversal symmetry
---------------------------------------------

All of the :math:`4 \times 4` matrices :math:`\boldsymbol{M}` described in Section `Syntax for the specification of one-electron operators`_ can be written as the Kronecker product of two :math:`2 \times 2` matrices :math:`\boldsymbol{A}` and :math:`\boldsymbol{B}`, such that :math:`\boldsymbol{M} = \boldsymbol{A} \otimes \boldsymbol{B}`.

In particular, :math:`\boldsymbol{A}` is one of the following real matrices: :math:`I_{2 \times 2}`, :math:`\sigma_x`, :math:`i \sigma_y`, and :math:`\sigma_z`, whereas :math:`\boldsymbol{B}` can be equal to :math:`I_{2 \times 2}`, :math:`i \sigma_x`, :math:`i \sigma_y`, or :math:`i \sigma_z`. Here, :math:`\sigma_x`, :math:`\sigma_y`, and :math:`\sigma_z` are the Pauli matrices, and :math:`i=\sqrt{-1}` is the imaginary unit. The individual matrix elements of :math:`\boldsymbol{M}` are given as

.. math::

    M_{pq} = A_{ij} \, B_{kl} \,,
    
where :math:`p=2(i-1)+k` and :math:`q=2(j-1)+l`.

In the DIRAC code, the structure of the Dirac equation with respect to time-reversal symmetry is displayed by using reordered Dirac 4-spinors, where the components are grouped on spin labels (:math:`\alpha`, :math:`\beta`) rather than large (upper) and small (lower) components (:math:`L`, :math:`S`), according to Eq. (7) of Ref. :cite:`Saue1999`. The same reordering applies to all the matrices :math:`\boldsymbol{M}`, which are then transformed accordingly.

In the original ordering, these matrices are given as

.. math::

    \boldsymbol{M} = \boldsymbol{A} \otimes \boldsymbol{B} =
    \begin{pmatrix}
    M_{L\alpha,L\alpha} & M_{L\alpha,L\beta} & M_{L\alpha,S\alpha} & M_{L\alpha,S\beta} \\
    M_{L\beta,L\alpha}  & M_{L\beta,L\beta}  & M_{L\beta,S\alpha}  & M_{L\beta,S\beta}  \\
    M_{S\alpha,L\alpha} & M_{S\alpha,L\beta} & M_{S\alpha,S\alpha} & M_{S\alpha,S\beta} \\
    M_{S\beta,L\alpha}  & M_{S\beta,L\beta}  & M_{S\beta,S\alpha}  & M_{S\beta,S\beta}
    \end{pmatrix}.

However, after reordering one clearly has the transformed matrices

.. math::

    \boldsymbol{M}^\prime = \boldsymbol{B} \otimes \boldsymbol{A} = 
    \begin{pmatrix}
    M_{L\alpha,L\alpha} & M_{L\alpha,S\alpha} & M_{L\alpha,L\beta} & M_{L\alpha,S\beta} \\
    M_{S\alpha,L\alpha} & M_{S\alpha,S\alpha} & M_{S\alpha,L\beta} & M_{S\alpha,S\beta} \\
    M_{L\beta,L\alpha}  & M_{L\beta,S\alpha}  & M_{L\beta,L\beta}  & M_{L\beta,S\beta}  \\
    M_{S\beta,L\alpha}  & M_{S\beta,S\alpha}  & M_{S\beta,L\beta}  & M_{S\beta,S\beta}
    \end{pmatrix}.

In other words, :math:`M^\prime_{PQ} = B_{kl} \, A_{ij}`, with :math:`P=2(k-1)+i` and :math:`Q=2(l-1)+j`. This reordering guaranties that all the matrices :math:`\boldsymbol{M}^\prime` are time-reversal symmetric.

Additionally, following the arguments given in Ref. :cite:`Saue1999` it is possible to make the following mapping, that involves the quaternion units: :math:`I_{2 \times 2} \rightarrow 1`, :math:`i \sigma_z \rightarrow \check{\imath}`, :math:`i \sigma_y \rightarrow \check{\jmath}`, and :math:`i \sigma_x \rightarrow \check{k}`. In this way, the matrices :math:`\boldsymbol{M}^\prime` can be mapped to :math:`2 \times 2` matrices that are linear combinations of the quaternion units.

The relationship between :math:`\boldsymbol{M}` and the quaternion form of the :math:`2 \times 2` mapped matrices (arising from :math:`\boldsymbol{M}^\prime`) is given in the following Table:

+------------------------+-----------------------------------------------------+-----------------------------------------------------+
| :math:`\boldsymbol{M}` | :math:`\boldsymbol{A} \otimes \boldsymbol{B}`       | Quaternion form of :math:`\boldsymbol{M}^\prime`    |
+========================+=====================================================+=====================================================+
| :math:`I_{4 \times 4}` | :math:`I_{2 \times 2} \otimes I_{2 \times 2}`       | :math:`I_{2 \times 2}`                              |
+------------------------+-----------------------------------------------------+-----------------------------------------------------+
| :math:`i\Sigma_z`      | :math:`I_{2 \times 2} \otimes i \sigma_z`           | :math:`\check{\imath}\, I_{2 \times 2}`             |
+------------------------+-----------------------------------------------------+-----------------------------------------------------+
| :math:`i\Sigma_y`      | :math:`I_{2 \times 2} \otimes i \sigma_y`           | :math:`\check{\jmath}\, I_{2 \times 2}`             |
+------------------------+-----------------------------------------------------+-----------------------------------------------------+
| :math:`i\Sigma_x`      | :math:`I_{2 \times 2} \otimes i \sigma_x`           | :math:`\check{k}\, I_{2 \times 2}`                  |
+------------------------+-----------------------------------------------------+-----------------------------------------------------+
| :math:`\beta`          | :math:`\sigma_z \otimes I_{2 \times 2}`             | :math:`\sigma_z`                                    |
+------------------------+-----------------------------------------------------+-----------------------------------------------------+
| :math:`i\beta\Sigma_z` | :math:`\sigma_z \otimes i \sigma_z`                 | :math:`\check{\imath}\, \sigma_z`                   |
+------------------------+-----------------------------------------------------+-----------------------------------------------------+
| :math:`i\beta\Sigma_y` | :math:`\sigma_z \otimes i \sigma_y`                 | :math:`\check{\jmath}\, \sigma_z`                   |
+------------------------+-----------------------------------------------------+-----------------------------------------------------+
| :math:`i\beta\Sigma_x` | :math:`\sigma_z \otimes i \sigma_x`                 | :math:`\check{k}\, \sigma_z`                        |
+------------------------+-----------------------------------------------------+-----------------------------------------------------+
| :math:`\beta\gamma^5`  | :math:`(i\sigma_y) \otimes I_{2 \times 2}`          | :math:`(i\sigma_y)`                                 |
+------------------------+-----------------------------------------------------+-----------------------------------------------------+
| :math:`i\beta\alpha_z` | :math:`(i\sigma_y) \otimes i \sigma_z`              | :math:`\check{\imath}\, (i\sigma_y)`                |
+------------------------+-----------------------------------------------------+-----------------------------------------------------+
| :math:`i\beta\alpha_y` | :math:`(i\sigma_y) \otimes i \sigma_y`              | :math:`\check{\jmath}\, (i\sigma_y)`                |
+------------------------+-----------------------------------------------------+-----------------------------------------------------+
| :math:`i\beta\alpha_x` | :math:`(i\sigma_y) \otimes i \sigma_x`              | :math:`\check{k}\, (i\sigma_y)`                     |
+------------------------+-----------------------------------------------------+-----------------------------------------------------+
| :math:`\gamma^5`       | :math:`\sigma_x \otimes I_{2 \times 2}`             | :math:`\sigma_x`                                    |
+------------------------+-----------------------------------------------------+-----------------------------------------------------+
| :math:`i\alpha_z`      | :math:`\sigma_x \otimes i \sigma_z`                 | :math:`\check{\imath}\, \sigma_x`                   |
+------------------------+-----------------------------------------------------+-----------------------------------------------------+
| :math:`i\alpha_y`      | :math:`\sigma_x \otimes i \sigma_y`                 | :math:`\check{\jmath}\, \sigma_x`                   |
+------------------------+-----------------------------------------------------+-----------------------------------------------------+
| :math:`i\alpha_x`      | :math:`\sigma_x \otimes i \sigma_x`                 | :math:`\check{k}\, \sigma_x`                        |
+------------------------+-----------------------------------------------------+-----------------------------------------------------+

It is important to note that although the time-reversal antisymmetric (T-odd) operators :math:`\boldsymbol{\alpha}`, :math:`\boldsymbol{\Sigma}`, :math:`\boldsymbol{\beta\alpha}`, and :math:`\boldsymbol{\beta\Sigma}` become time-reversal symmetric (T-even) once they are reordered and written in quaternion form as described above, DIRAC saves the information that these operators are T-odd in the original ordering.

For more details about the quaternion symmetry scheme used in DIRAC, please refer to :cite:`Saue1999`.

Note: For complete consistency, perhaps the name of the operators :math:`i\boldsymbol{\alpha}`, :math:`i\boldsymbol{\Sigma}`, :math:`i\boldsymbol{\beta\alpha}`, and :math:`i\boldsymbol{\beta\Sigma}` (see `Operator type keywords`_) should have been ``XiALPHA``, ``YiALPHA``, ``ZiALPHA``, ``XiAVECTOR``, ``YiAVECTOR``, ``ZiAVECTOR``, ``iALPHADOT``, ``XiBETAALP``, ``YiBETAALP``, ``ZiBETAALP``, ``XiSIGMA``, ``YiSIGMA``, ``ZiSIGMA``, ``XiBETASIG``, ``YiBETASIG``, ``ZiBETASIG``, and ``iSIGMADOT``. For historical reasons, we do not use this convention.

Note for advanced users: The operator iBETAGAMMA5 is implemented in the code just as :math:`\beta\gamma^5`, but its time-reversal symmetry without imaginary phase is assigned to be odd.
