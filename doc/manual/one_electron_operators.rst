:orphan:


.. _one_electron_operators:

One-electron operators
======================

Syntax for specifying one-electron operators
--------------------------------------------

In DIRAC, a general one-electron 4-component operator is generated from linear combinations of operators having the basic form

.. math::

  \hat{\boldsymbol{P}} = f \, \boldsymbol{M} \, \hat{\Omega},

where :math:`f` is a scalar factor, :math:`\hat{\Omega}` is a scalar operator, and :math:`\boldsymbol{M}` is one of the following :math:`4 \times 4` time-reversal symmetric matrices:

.. math::

  & \boldsymbol{I}, \, \boldsymbol{\beta}, \, \boldsymbol{\gamma}_5, \, \boldsymbol{\beta}\boldsymbol{\gamma}_5,

  & i\boldsymbol{\alpha}_x, \, i\boldsymbol{\alpha}_y, \, i\boldsymbol{\alpha}_z, \, i\boldsymbol{\beta}\boldsymbol{\alpha}_x, \, i\boldsymbol{\beta}\boldsymbol{\alpha}_y, \, i\boldsymbol{\beta}\boldsymbol{\alpha}_z,

  & i\boldsymbol{\Sigma}_x, \, i\boldsymbol{\Sigma}_y, \, i\boldsymbol{\Sigma}_z, \, i\boldsymbol{\beta}\boldsymbol{\Sigma}_x, \, i\boldsymbol{\beta}\boldsymbol{\Sigma}_y, \, i\boldsymbol{\beta}\boldsymbol{\Sigma}_z.

One thing to notice is that :math:`\boldsymbol{I}`, :math:`\boldsymbol{\beta}`, :math:`\boldsymbol{\gamma}_5`, and :math:`\boldsymbol{\beta}\boldsymbol{\gamma}_5` are time-reversal symmetric operators, and they are implemented with such symmetry in DIRAC. On the other hand, :math:`\boldsymbol{\alpha}_k`, :math:`\boldsymbol{\Sigma}_k`, :math:`\boldsymbol{\beta\alpha}_k`, and :math:`\boldsymbol{\beta\Sigma}_k` (with :math:`k=x,y`, and :math:`z`) are time-reversal antisymmetric operators. This information is saved by the program, but these matrices (i.e., :math:`\boldsymbol{\alpha}_k`, :math:`\boldsymbol{\Sigma}_k`, :math:`\boldsymbol{\beta\alpha}_k`, and :math:`\boldsymbol{\beta\Sigma}_k`) are multiplied by the imaginary unit :math:`i=\sqrt{-1}` in order to make them time-reversal symmetric and hence fit into the quaternion symmetry scheme of DIRAC (see :cite:`Saue1999` and :cite:`Salek2005` for more information).

Another important aspect to note is that the :math:`4 \times 4` matrices :math:`\boldsymbol{M}` used in DIRAC correspond to reordered Dirac 4-spinors, as described in Section `Quaternion algebra and time-reversal symmetry`_.

Further details about the implemented operators are provided in Section `Quaternion algebra and time-reversal symmetry`_, where the connection between the Dirac matrices and the quaternion algebra is given.


Operator specification
^^^^^^^^^^^^^^^^^^^^^^

Operators are specified by the keyword :ref:`HAMILTONIAN_.OPERATOR` with the following arguments::

  .OPERATOR
   'user defined operator name'
   operator type keyword
   one-electron scalar operator label for each component
   CMULT
   FACTORS
   factor for each component
   COMFACTOR
   common factor for all components

Note that the arguments following the keyword :ref:`HAMILTONIAN_.OPERATOR` must start with a blank.

The arguments are optional, except for the one-electron scalar operator label. If no operator type keyword is given, DIRAC assumes the operator type to be :math:`\boldsymbol{I}`, the diagonal :math:`4 \times 4` matrix. Component factors under the keyword ``FACTORS``, as well as the common factor under the keyword ``COMFACTOR``, are all equal to one if not specified.

The keyword ``CMULT`` assures multiplication of the individual factors (or common factor) by the speed of light in vacuum, :math:`c`. This option has the further advantage that ``CMULT`` follows any user-specified modification of the speed of light, as provided by :ref:`GENERAL_.CVALUE`.


.. _Operator type keywords:

Operator type keywords
^^^^^^^^^^^^^^^^^^^^^^

There are 22 basic operator types implemented in DIRAC. They are listed in the following Table, containing their keywords, operator forms, associated number of factors (NF), and time-reversal symmetry (TRS) without imaginary phase:

===========  =========================================================================================================================================  ======  ===================
**Keyword**  **Operator form**                                                                                                                          **NF**  **TRS (without i)**
===========  =========================================================================================================================================  ======  ===================
DIAGONAL     :math:`f \boldsymbol{I} \hat{\Omega}`                                                                                                      1       Symmetric
BETA         :math:`f \boldsymbol{\beta} \hat{\Omega}`                                                                                                  1       Symmetric
GAMMA5       :math:`f \boldsymbol{\gamma}_5 \hat{\Omega}`                                                                                               1       Symmetric
BETAGAMM     :math:`f \boldsymbol{\beta} \boldsymbol{\gamma}_5 \hat{\Omega}`                                                                            1       Symmetric
XALPHA       :math:`f i\boldsymbol{\alpha}_x \hat{\Omega}`                                                                                              1       Antisymmetric
YALPHA       :math:`f i\boldsymbol{\alpha}_y \hat{\Omega}`                                                                                              1       Antisymmetric
ZALPHA       :math:`f i\boldsymbol{\alpha}_z \hat{\Omega}`                                                                                              1       Antisymmetric
XAVECTOR     :math:`f_1 i\boldsymbol{\alpha}_y \hat{\Omega}_z - f_2 i\boldsymbol{\alpha}_z \hat{\Omega}_y`                                              2       Antisymmetric
YAVECTOR     :math:`f_1 i\boldsymbol{\alpha}_z \hat{\Omega}_x - f_2 i\boldsymbol{\alpha}_x \hat{\Omega}_z`                                              2       Antisymmetric
ZAVECTOR     :math:`f_1 i\boldsymbol{\alpha}_x \hat{\Omega}_y - f_2 i\boldsymbol{\alpha}_y \hat{\Omega}_x`                                              2       Antisymmetric
ALPHADOT     :math:`f_1 i\boldsymbol{\alpha}_x \hat{\Omega}_x + f_2 i\boldsymbol{\alpha}_y \hat{\Omega}_y + f_3 i\boldsymbol{\alpha}_z \hat{\Omega}_z`  3       Antisymmetric
XBETAALP     :math:`f i\boldsymbol{\beta} \boldsymbol{\alpha}_x \hat{\Omega}`                                                                           1       Antisymmetric
YBETAALP     :math:`f i\boldsymbol{\beta} \boldsymbol{\alpha}_y \hat{\Omega}`                                                                           1       Antisymmetric
ZBETAALP     :math:`f i\boldsymbol{\beta} \boldsymbol{\alpha}_z \hat{\Omega}`                                                                           1       Antisymmetric
XSIGMA       :math:`f i\boldsymbol{\Sigma}_x \hat{\Omega}`                                                                                              1       Antisymmetric
YSIGMA       :math:`f i\boldsymbol{\Sigma}_y \hat{\Omega}`                                                                                              1       Antisymmetric
ZSIGMA       :math:`f i\boldsymbol{\Sigma}_z \hat{\Omega}`                                                                                              1       Antisymmetric
SIGMADOT     :math:`f_1 i\boldsymbol{\Sigma}_x \hat{\Omega}_x + f_2 i\boldsymbol{\Sigma}_y \hat{\Omega}_y + f_3 i\boldsymbol{\Sigma}_z \hat{\Omega}_z`  1       Antisymmetric
XBETASIG     :math:`f i\boldsymbol{\beta} \boldsymbol{\Sigma}_x \hat{\Omega}`                                                                           1       Antisymmetric
YBETASIG     :math:`f i\boldsymbol{\beta} \boldsymbol{\Sigma}_y \hat{\Omega}`                                                                           1       Antisymmetric
ZBETASIG     :math:`f i\boldsymbol{\beta} \boldsymbol{\Sigma}_z \hat{\Omega}`                                                                           1       Antisymmetric
iBETAGAM     :math:`f \boldsymbol{\beta} \boldsymbol{\gamma}_5 \hat{\Omega}`                                                                            1       Antisymmetric
===========  =========================================================================================================================================  ======  ===================

`Note for advanced users:` iBETAGAM corresponds to the time-reversal symmetric operator :math:`f \boldsymbol{\beta} \boldsymbol{\gamma}_5 \hat{\Omega}`, but has been assigned a time-reversal antisymmetry. This operator can be used, for example, for calculations of the molecular enhancement factors of the electron electric dipole moment (eEDM). For further details, see Section `Quaternion algebra and time-reversal symmetry`_.


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
| FERMI C     | One-electron Fermi contact       | FC NAMab       | :math:`\hat{\Omega}_1 = \frac{4 \pi g_{e}}{3} \delta(\boldsymbol{r}-\boldsymbol{R}_K)` |
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

where :math:`\boldsymbol{R}_{CM}` and :math:`\boldsymbol{R}_K` are the position vectors of the nuclear center of mass and the nucleus :math:`K` with respect to the origin of coordinates, respectively. In addition, the operators :math:`x`, :math:`y`, and :math:`z` are the Cartesian components of the electron position vector operator with respect to the origin of coordinates.

Note that in the one-electron Fermi-contact integrals :math:`g_e` is the electronic g-factor (i.e., :math:`g_e = 2` unless QED effects are taken into account), 'NAM' are the three first letters in the name of atom :math:`K` as given in the .mol file, and 'ab' stands for the number of the symmetry-adapted nucleus.


Other one-electron operator labels are:

=========== ====================================================================================================================
**Keyword** **Description**
=========== ====================================================================================================================
ANGLON      Angular momentum around the nuclei
CARMOM      Cartesian moments (symmetric) :math:`(l + 1)(l + 2)/2` components (:math:`l = i + j + k`)
CM1         First order magnetic field derivatives of electric field
CM2         Second order magnetic field derivatives of electric field
DARWIN      Darwin-type integrals
DIASUS      Angular London orbital contribution to diamagnetic susceptibility
DSO         Diamagnetic spin-orbit integrals
DSUSCGO     Diamagnetic susceptibility with common gauge origin
DSUSLH      Angular London orbital contribution to diamagnetic susceptibility
DSUSNOL     Diamagnetic susceptibility without London contribution
ELFGRDC     Electric field gradient at the individual nuclei, cartesian
ELFGRDS     Electric field gradient at the individual nuclei, spherical
EXPIKR      Cosine and sine
HBDO        Half B-differentiated overlap matrix
HDO         Half-derivative overlap integrals
HDOBR       Ket-differentiation of HDO-integrals with respect to magnetic field
LONMOM      London orbital contribution to angular momentum
MAGMOM      One-electron contributions to magnetic moment
MASSVEL     Mass velocity
NEFIELD     Electric field at the individual nuclei
NSTCGO      Nuclear shielding with common gauge origin
NUCPOT      Potential energy of the interaction of electrons with individual nuclei, divided by the nuclear charge
NUCSHI      Nuclear shielding tensor
NUCSLO      London orbital contribution to nuclear shielding tensor
NUCSNLO     Nuclear shielding integrals without London orbital contribution
PSO         Paramagnetic spin-orbit
S1MAG       Second order contribution from overlap matrix to magnetic properties
S1MAGL      Bra-differentiation of overlap matrix with respect to magnetic field
S1MAGR      Ket-differentiation of overlap matrix with respect to magnetic field
SDFC        Spin-dipole + Fermi contact
SOLVENT     Electronic solvent integrals
SPHMOM      Spherical moments (real combinations), symmetric, (:math:`2l + 1`) components (:math:`m = 0, -1, +1, ..., -l, +l`)
SPIN-DI     Spin-dipole integrals
SPNORB      Spatial spin-orbit
SQHDO       Half-derivative overlap integrals not to be antisymmetrized
SQHDOR      Half-derivative overlap integrals not to be anti-symmetrized
SQOVLAP     Second order derivatives overlap integrals
=========== ====================================================================================================================


Examples of using various operators
-----------------------------------

We give here several concrete examples on how to construct operators for various properties.


Kinetic part of the Dirac Hamiltonian
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The kinetic part of the electronic Dirac Hamiltonian is given (in a.u.) by

.. math::

  \hat{H}^K & = c\vec{\boldsymbol{\alpha}}\cdot\hat{\vec{p}} = c\left(\boldsymbol{\alpha}_x\hat{p}_x+\boldsymbol{\alpha}_y\hat{p}_y+\boldsymbol{\alpha}_z\hat{p}_z\right)

            & = c \left[ i\boldsymbol{\alpha}_x\left(-\frac{\partial}{\partial x}\right) + i\boldsymbol{\alpha}_y\left(-\frac{\partial}{\partial y}\right) + i\boldsymbol{\alpha}_z\left(-\frac{\partial}{\partial z}\right) \right].

This operator may be specified by::

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


Magnetic interactions
^^^^^^^^^^^^^^^^^^^^^

Another example is the operator that describes the electromagnetic interaction between electrons and a uniform external magnetic field: :math:`\hat{H}^B=-\frac{c}{2}\vec{B}\cdot\vec{\boldsymbol{\alpha}}\times\hat{\vec{r}}`.

The operator to use in this case is :math:`-\frac{c}{2}\vec{\boldsymbol{\alpha}}\times\vec{r}`, and in DIRAC each component is calculated individually. For example, the :math:`x`-component of this operator is :math:`-\frac{c}{2}\left(\boldsymbol{\alpha}_y z - \boldsymbol{\alpha}_z y\right)`, and can be requested using::

  .OPERATOR
   'B_x'
   XAVECTOR
   ZDIPLEN
   YDIPLEN
   CMULT
   COMFACTOR
   -0.5

The program will assume all operators to be Hermitian. The operator ``XAVECTOR`` is defined as :math:`f_1 i\boldsymbol{\alpha}_y \hat{\Omega}_z-f_2 i\boldsymbol{\alpha}_z \hat{\Omega}_y` (see `Operator type keywords`_), and in this particular case the scalar operators are :math:`\hat{\Omega}_z=z` and :math:`\hat{\Omega}_y=y` (see ``ZDIPLEN`` and ``YDIPLEN`` in `One-electron scalar operator labels`_). Additionally, :math:`f_1=f_2=-\frac{c}{2}`. Then, DIRAC will assume the operator to be :math:`-\frac{c}{2}\left(i\boldsymbol{\alpha}_y z - i\boldsymbol{\alpha}_z y\right)`.

The imaginary phase :math:`i` makes the operator time-reversal symmetric.


Diagonal operators
^^^^^^^^^^^^^^^^^^

If no other argument is given, the program assumes the operator to be diagonal and expects the operator name to be the component label, for instance::

  .OPERATOR
   OVERLAP


Finite-field perturbation: Electric dipole moment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Another example is a finite-perturbation calculation where the perturbed Hamiltonian operator arising from the interaction between the electrons and a uniform external electric field oriented along the :math:`z`-axis (:math:`\hat{H}'=-q E_z \hat{z}=e E_z \hat{z}`) is added to the unperturbed Hamiltonian :math:`\hat{H}_0` defined in :ref:`**HAMILTONIAN` (default: Dirac-Coulomb).

If we take :math:`e E_z = 0.01` (in a.u.), the :math:`\hat{z}` dipole length operator will be used to build the total Hamiltonian :math:`\hat{H} = \hat{H}_0 + 0.01 \, \hat{z}` by adding the perturbed Hamiltonian :math:`\hat{H}'` to the Dirac-Coulomb Hamiltonian :math:`\hat{H}_0` using::

  **HAMILTONIAN
  .OPERATOR
   ZDIPLEN
   COMFACTOR
   0.01

`Note:` Don't forget to decrease the symmetry of your system.


Fermi-contact integrals
^^^^^^^^^^^^^^^^^^^^^^^

Here is an example where the Fermi-contact (FC) integrals for a certain nucleus :math:`K`
are added to the Hamiltonian in a finite-field calculation.

Let's assume we are looking at a PbX dimer (order in the .mol file: 1. Pb, 2. X) and we want to add (under :ref:`**HAMILTONIAN`) the FC interaction for the Pb nucleus to the unperturbed Hamiltonian :math:`\hat{H}_0` as a perturbation with a given field strength of :math:`10^{-9}` a.u..

The FC integrals are called using the label ``'FC NAMab'`` (see `One-electron scalar operator labels`_), and the total perturbed Hamiltonian :math:`\hat{H}=\hat{H}_0+\hat{H}'` (where, in a.u., :math:`\hat{H}'=10^{-9} \, \frac{4\pi}{3} \, g_e \delta(\boldsymbol{r}-\boldsymbol{R}_\text{Pb})`) can be requested using::

  **HAMILTONIAN
  .OPERATOR
   'Density at Pb nucleus'
   DIAGONAL
   'FC Pb 01'
   FACTORS
   0.000000001

**Important note:** The values obtained after the fit of the finite-field energies need to be scaled by :math:`\frac{3}{4 \pi g_{e}} = \frac{1}{8.3872954891254192}` in order to get the raw density :math:`\delta(\boldsymbol{r}-\boldsymbol{R}_\text{Pb})`.


The next example shows how to calculate the electronic density at the different nuclei in the molecule, as the expectation values :math:`\langle -\delta(\boldsymbol{r}-\boldsymbol{R}_K) \rangle` (for each nucleus :math:`K` in the system) for a Dirac-Coulomb Hartree-Fock wave function including a decomposition of the molecular orbital contributions to the density::

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

For further details, see :ref:`PROPERTIES_.RHONUC`.


Cartesian moment expectation value
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In the following example we calculate the cartesian multipole moment expectation value :math:`\langle x^1 y^2 z^3 \rangle` for a Levy-Leblond Hartree-Fock wave function::

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

In particular, :math:`\boldsymbol{A}` is one of the real matrices :math:`I_{2 \times 2}`, :math:`\sigma_x`, :math:`i \sigma_y`, or :math:`\sigma_z`, whereas :math:`\boldsymbol{B}` is :math:`I_{2 \times 2}`, :math:`i \sigma_x`, :math:`i \sigma_y`, or :math:`i \sigma_z`. Here, :math:`\sigma_x`, :math:`\sigma_y`, and :math:`\sigma_z` are the Pauli matrices, and :math:`i` is the imaginary unit. The individual matrix elements of :math:`\boldsymbol{M}` are given as

.. math::

    M_{pq} = A_{ij} \, B_{kl} \, ,

where :math:`p=2(i-1)+k` and :math:`q=2(j-1)+l`.

In the DIRAC code, the structure of the Dirac equation with respect to time-reversal symmetry is displayed by using reordered Dirac 4-spinors, where the components are grouped on spin labels (:math:`\alpha`, :math:`\beta`) rather than large (upper) and small (lower) components (:math:`L`, :math:`S`),

.. math::

   \begin{pmatrix}
   \psi^{L\alpha} \\ \psi^{L\beta} \\ \psi^{S\alpha} \\ \psi^{S\beta}
   \end{pmatrix}
   \rightarrow
   \begin{pmatrix}
   \psi^{L\alpha} \\ \psi^{S\alpha} \\ \psi^{L\beta} \\ \psi^{S\beta}
   \end{pmatrix}.

The same reordering applies to all the matrices :math:`\boldsymbol{M}`, which are then transformed accordingly. In the original ordering, these matrices are given as

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

In other words, :math:`M^\prime_{PQ} = B_{kl} \, A_{ij}`, with :math:`P=2(k-1)+i` and :math:`Q=2(l-1)+j`. This reordering guarantees that all the matrices :math:`\boldsymbol{M}^\prime` are time-reversal symmetric.

Additionally, following the arguments given in Ref. :cite:`Saue1999`, it is possible to make the following mapping involving the quaternion units: :math:`I_{2 \times 2} \rightarrow 1`, :math:`i \sigma_z \rightarrow \check{\imath}`, :math:`i \sigma_y \rightarrow \check{\jmath}`, and :math:`i \sigma_x \rightarrow \check{k}`. Thus, the matrices :math:`\boldsymbol{M}^\prime` can be represented as :math:`2 \times 2` matrices, which are linear combinations of the quaternion units.

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

For further details on the quaternion symmetry scheme used in DIRAC, see Ref. :cite:`Saue1999`.

Note: For complete consistency, the names of the operators could have been :math:`i\boldsymbol{\alpha}`, :math:`i\boldsymbol{\Sigma}`, :math:`i\boldsymbol{\beta\alpha}`, and :math:`i\boldsymbol{\beta\Sigma}` (see `Operator type keywords`_) should have been ``XiALPHA``, ``YiALPHA``, ``ZiALPHA``, ``XiAVECTOR``, ``YiAVECTOR``, ``ZiAVECTOR``, ``iALPHADOT``, ``XiBETAALP``, ``YiBETAALP``, ``ZiBETAALP``, ``XiSIGMA``, ``YiSIGMA``, ``ZiSIGMA``, ``XiBETASIG``, ``YiBETASIG``, ``ZiBETASIG``, and ``iSIGMADOT``. However, due to historical reasons, this convention is not used in DIRAC.

Note for advanced users: The operator iBETAGAM is implemented in the code just as :math:`\beta\gamma^5`, but its time-reversal symmetry without imaginary phase is designated as odd.
