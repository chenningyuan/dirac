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

  & \boldsymbol{I}, \, \boldsymbol{\beta}, \, \boldsymbol{\gamma}^5, \, \boldsymbol{\beta}\boldsymbol{\gamma}^5,

  & i\boldsymbol{\alpha}_x, \, i\boldsymbol{\alpha}_y, \, i\boldsymbol{\alpha}_z, \, i\boldsymbol{\beta}\boldsymbol{\alpha}_x, \, i\boldsymbol{\beta}\boldsymbol{\alpha}_y, \, i\boldsymbol{\beta}\boldsymbol{\alpha}_z,

  & i\boldsymbol{\Sigma}_x, \, i\boldsymbol{\Sigma}_y, \, i\boldsymbol{\Sigma}_z, \, i\boldsymbol{\beta}\boldsymbol{\Sigma}_x, \, i\boldsymbol{\beta}\boldsymbol{\Sigma}_y, \, i\boldsymbol{\beta}\boldsymbol{\Sigma}_z.

One thing to notice is that :math:`\boldsymbol{I}`, :math:`\boldsymbol{\beta}`, :math:`\boldsymbol{\gamma}^5`, and :math:`\boldsymbol{\beta}\boldsymbol{\gamma}^5` are time-reversal symmetric operators, and they are implemented with such symmetry in DIRAC. On the other hand, :math:`\boldsymbol{\alpha}_k`, :math:`\boldsymbol{\Sigma}_k`, :math:`\boldsymbol{\beta\alpha}_k`, and :math:`\boldsymbol{\beta\Sigma}_k` (with :math:`k=x,y`, and :math:`z`) are time-reversal antisymmetric operators. This information is saved by the program, but these matrices (i.e., :math:`\boldsymbol{\alpha}_k`, :math:`\boldsymbol{\Sigma}_k`, :math:`\boldsymbol{\beta\alpha}_k`, and :math:`\boldsymbol{\beta\Sigma}_k`) are multiplied by the imaginary unit :math:`i=\sqrt{-1}` in order to make them time-reversal symmetric and hence fit into the quaternion symmetry scheme of DIRAC (see :cite:`Saue1999` and :cite:`Salek2005` for more information).

Another important aspect to note is that DIRAC does not use the :math:`4 \times 4` :math:`\boldsymbol{M}` matrices in the standard representation. These matrices are reordered following a criterion described in Section `Quaternion algebra and time-reversal symmetry`_. Further details about the implemented operators are provided in that Section, where the connection between the Dirac matrices and the quaternion algebra is given.


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

There are 23 basic operator types implemented in DIRAC. They are listed in the following Table, containing their keywords, operator forms, associated number of factors (NF), and time-reversal symmetry (TRS) without imaginary phase:

===========  ===============================================================================================================================================================================================  ======  ===================
**Keyword**  **Operator form**                                                                                                                                                                                **NF**  **TRS (without i)**
===========  ===============================================================================================================================================================================================  ======  ===================
DIAGONAL     :math:`f \boldsymbol{I} \hat{\Omega}`                                                                                                                                                            1       Symmetric
BETA         :math:`f \boldsymbol{\beta} \hat{\Omega}`                                                                                                                                                        1       Symmetric
GAMMA5       :math:`f \boldsymbol{\gamma}^5 \hat{\Omega}`                                                                                                                                                     1       Symmetric
BETAGAMM     :math:`f \boldsymbol{\beta} \boldsymbol{\gamma}^5 \hat{\Omega}`                                                                                                                                  1       Symmetric
XALPHA       :math:`f i\boldsymbol{\alpha}_x \hat{\Omega}`                                                                                                                                                    1       Antisymmetric
YALPHA       :math:`f i\boldsymbol{\alpha}_y \hat{\Omega}`                                                                                                                                                    1       Antisymmetric
ZALPHA       :math:`f i\boldsymbol{\alpha}_z \hat{\Omega}`                                                                                                                                                    1       Antisymmetric
XAVECTOR     :math:`f_1 i\boldsymbol{\alpha}_y \hat{\Omega}_z - f_2 i\boldsymbol{\alpha}_z \hat{\Omega}_y`                                                                                                    2       Antisymmetric
YAVECTOR     :math:`f_1 i\boldsymbol{\alpha}_z \hat{\Omega}_x - f_2 i\boldsymbol{\alpha}_x \hat{\Omega}_z`                                                                                                    2       Antisymmetric
ZAVECTOR     :math:`f_1 i\boldsymbol{\alpha}_x \hat{\Omega}_y - f_2 i\boldsymbol{\alpha}_y \hat{\Omega}_x`                                                                                                    2       Antisymmetric
ALPHADOT     :math:`f_1 i\boldsymbol{\alpha}_x \hat{\Omega}_x + f_2 i\boldsymbol{\alpha}_y \hat{\Omega}_y + f_3 i\boldsymbol{\alpha}_z \hat{\Omega}_z`                                                        3       Antisymmetric
XBETAALP     :math:`f i\boldsymbol{\beta} \boldsymbol{\alpha}_x \hat{\Omega}`                                                                                                                                 1       Antisymmetric
YBETAALP     :math:`f i\boldsymbol{\beta} \boldsymbol{\alpha}_y \hat{\Omega}`                                                                                                                                 1       Antisymmetric
ZBETAALP     :math:`f i\boldsymbol{\beta} \boldsymbol{\alpha}_z \hat{\Omega}`                                                                                                                                 1       Antisymmetric
XSIGMA       :math:`f i\boldsymbol{\Sigma}_x \hat{\Omega}`                                                                                                                                                    1       Antisymmetric
YSIGMA       :math:`f i\boldsymbol{\Sigma}_y \hat{\Omega}`                                                                                                                                                    1       Antisymmetric
ZSIGMA       :math:`f i\boldsymbol{\Sigma}_z \hat{\Omega}`                                                                                                                                                    1       Antisymmetric
SIGMADOT     :math:`f_1 i\boldsymbol{\Sigma}_x \hat{\Omega}_x + f_2 i\boldsymbol{\Sigma}_y \hat{\Omega}_y + f_3 i\boldsymbol{\Sigma}_z \hat{\Omega}_z`                                                        1       Antisymmetric
XBETASIG     :math:`f i\boldsymbol{\beta} \boldsymbol{\Sigma}_x \hat{\Omega}`                                                                                                                                 1       Antisymmetric
YBETASIG     :math:`f i\boldsymbol{\beta} \boldsymbol{\Sigma}_y \hat{\Omega}`                                                                                                                                 1       Antisymmetric
ZBETASIG     :math:`f i\boldsymbol{\beta} \boldsymbol{\Sigma}_z \hat{\Omega}`                                                                                                                                 1       Antisymmetric
BEALDOT      :math:`f_1 i\boldsymbol{\beta}\boldsymbol{\alpha}_x \hat{\Omega}_x + f_2 i\boldsymbol{\beta}\boldsymbol{\alpha}_y \hat{\Omega}_y + f_3 i\boldsymbol{\beta}\boldsymbol{\alpha}_z \hat{\Omega}_z`  3       Antisymmetric
iBETAGAM     :math:`f \boldsymbol{\beta} \boldsymbol{\gamma}^5 \hat{\Omega}`                                                                                                                                  1       Antisymmetric
===========  ===============================================================================================================================================================================================  ======  ===================

`Note for advanced users:` iBETAGAM corresponds to the time-reversal symmetric operator :math:`f \boldsymbol{\beta} \boldsymbol{\gamma}^5 \hat{\Omega}`, but has been assigned a time-reversal antisymmetry. This operator can be used, for example, for calculations of the molecular enhancement factors of the electron electric dipole moment (eEDM). For further details, see Section `Quaternion algebra and time-reversal symmetry`_.


.. _One-electron scalar operator labels:

One-electron scalar operator labels
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

+-------------+----------------------------------+-------------+----------------+------------------------------------------------------------------------------------------------+
| **Operator**| **Integral description**         | **Symmetry**| **Components** | **Operators**                                                                                  |
| **label**   |                                  |             |                |                                                                                                |
+=============+==================================+=============+================+================================================================================================+
| KINENER     | NR kinetic energy                |     S       | KINENER        | :math:`\hat{\Omega}_1 = -\frac{1}{2}\hat{\nabla}^2`                                            |
+-------------+----------------------------------+-------------+----------------+------------------------------------------------------------------------------------------------+
| MOLFIELD    | Nuclear attraction               |     S       | MOLFIELD       | :math:`\hat{\Omega}_1 = \sum_K Z_K |\boldsymbol{r}-\boldsymbol{R}_K|^{-1}`                     |
+-------------+----------------------------------+-------------+----------------+------------------------------------------------------------------------------------------------+
| OVERLAP     | Overlap                          |     S       | OVERLAP        | :math:`\hat{\Omega}_1 = 1`                                                                     |
+-------------+----------------------------------+-------------+----------------+------------------------------------------------------------------------------------------------+
| FERMI C     | Fermi contact                    |     S       | FC NAMab       | :math:`\hat{\Omega}_1 = \frac{4 \pi g_{e}}{3} \delta(\boldsymbol{r}-\boldsymbol{R}_K)`         |
+-------------+----------------------------------+-------------+----------------+------------------------------------------------------------------------------------------------+
| PVC         | Gaussian nuclear density         |     S       | PVCNAMab       | :math:`\hat{\Omega}_1 = \rho^G_K(\boldsymbol{r})`                                              |
|             | normalized to 1,                 |             |                |                :math:`= (\zeta_K/\pi)^{3/2} \exp(-\zeta_K |\boldsymbol{r}-\boldsymbol{R}_K|^2)`|
|             | :math:`\rho^G_K(\boldsymbol{r})` |             |                |                                                                                                |
+-------------+----------------------------------+-------------+----------------+------------------------------------------------------------------------------------------------+
| GRAD FC     | Gradient of the Dirac delta      |     S       | XGFNAMab       | :math:`\hat{\Omega}_1 = \frac{\partial}{\partial x} \delta(\boldsymbol{r}-\boldsymbol{R}_K)`   |
|             | distribution                     |             +----------------+------------------------------------------------------------------------------------------------+
|             |                                  |             | YGFNAMab       | :math:`\hat{\Omega}_1 = \frac{\partial}{\partial y} \delta(\boldsymbol{r}-\boldsymbol{R}_K)`   |
|             |                                  |             +----------------+------------------------------------------------------------------------------------------------+
|             |                                  |             | ZGFNAMab       | :math:`\hat{\Omega}_1 = \frac{\partial}{\partial z} \delta(\boldsymbol{r}-\boldsymbol{R}_K)`   |
+-------------+----------------------------------+-------------+----------------+------------------------------------------------------------------------------------------------+
| DRHO        | Gradient of the Gaussian nuclear |     S       | DRHO lmn       | :math:`\hat{\Omega}_1 = \frac{\partial}{\partial x} \rho^G_K(\boldsymbol{r})`                  |
|             | density                          |             +----------------+------------------------------------------------------------------------------------------------+
|             | :math:`\rho^G_K(\boldsymbol{r})` |             | DRHO opq       | :math:`\hat{\Omega}_1 = \frac{\partial}{\partial y} \rho^G_K(\boldsymbol{r})`                  |
|             |                                  |             +----------------+------------------------------------------------------------------------------------------------+
|             |                                  |             | DRHO rst       | :math:`\hat{\Omega}_1 = \frac{\partial}{\partial z} \rho^G_K(\boldsymbol{r})`                  |
+-------------+----------------------------------+-------------+----------------+------------------------------------------------------------------------------------------------+
| ANGMOM      | Orbital angular momentum         |     A       | XANGMOM        | :math:`\hat{\Omega}_1 = -[(\boldsymbol{r}-\boldsymbol{R}_{CM})\times\hat{\nabla}]_x`           |
|             | around CM                        |             +----------------+------------------------------------------------------------------------------------------------+
|             |                                  |             | YANGMOM        | :math:`\hat{\Omega}_2 = -[(\boldsymbol{r}-\boldsymbol{R}_{CM})\times\hat{\nabla}]_y`           |
|             |                                  |             +----------------+------------------------------------------------------------------------------------------------+
|             |                                  |             | ZANGMOM        | :math:`\hat{\Omega}_3 = -[(\boldsymbol{r}-\boldsymbol{R}_{CM})\times\hat{\nabla}]_z`           |
+-------------+----------------------------------+-------------+----------------+------------------------------------------------------------------------------------------------+
| DIPLEN      | Dipole length                    |     S       | XDIPLEN        | :math:`\hat{\Omega}_1 = x`                                                                     |
|             |                                  |             +----------------+------------------------------------------------------------------------------------------------+
|             |                                  |             | YDIPLEN        | :math:`\hat{\Omega}_2 = y`                                                                     |
|             |                                  |             +----------------+------------------------------------------------------------------------------------------------+
|             |                                  |             | ZDIPLEN        | :math:`\hat{\Omega}_3 = z`                                                                     |
+-------------+----------------------------------+-------------+----------------+------------------------------------------------------------------------------------------------+
| DIPVEL      | Dipole velocity                  |     A       | XDIPVEL        | :math:`\hat{\Omega}_1 = \frac{\partial}{\partial x}`                                           |
|             |                                  |             +----------------+------------------------------------------------------------------------------------------------+
|             |                                  |             | YDIPVEL        | :math:`\hat{\Omega}_2 = \frac{\partial}{\partial y}`                                           |
|             |                                  |             +----------------+------------------------------------------------------------------------------------------------+
|             |                                  |             | ZDIPVEL        | :math:`\hat{\Omega}_3 = \frac{\partial}{\partial z}`                                           |
+-------------+----------------------------------+-------------+----------------+------------------------------------------------------------------------------------------------+
| QUADRUP     | Quadrupole moment                |     S       | XXQUADRU       | :math:`\hat{\Omega}_1 = \frac{1}{4} (x^2-r^2)`                                                 |
|             |                                  |             +----------------+------------------------------------------------------------------------------------------------+
|             |                                  |             | XYQUADRU       | :math:`\hat{\Omega}_2 = \frac{1}{4} xy`                                                        |
|             |                                  |             +----------------+------------------------------------------------------------------------------------------------+
|             |                                  |             | XZQUADRU       | :math:`\hat{\Omega}_3 = \frac{1}{4} xz`                                                        |
|             |                                  |             +----------------+------------------------------------------------------------------------------------------------+
|             |                                  |             | YYQUADRU       | :math:`\hat{\Omega}_4 = \frac{1}{4} (y^2-r^2)`                                                 |
|             |                                  |             +----------------+------------------------------------------------------------------------------------------------+
|             |                                  |             | YZQUADRU       | :math:`\hat{\Omega}_5 = \frac{1}{4} yz`                                                        |
|             |                                  |             +----------------+------------------------------------------------------------------------------------------------+
|             |                                  |             | ZZQUADRU       | :math:`\hat{\Omega}_6 = \frac{1}{4} (z^2-r^2)`                                                 |
+-------------+----------------------------------+-------------+----------------+------------------------------------------------------------------------------------------------+
| SECMOM      | Second moment                    |     S       | XXSECMOM       | :math:`\hat{\Omega}_1 = - x^2`                                                                 |
|             |                                  |             +----------------+------------------------------------------------------------------------------------------------+
|             |                                  |             | XYSECMOM       | :math:`\hat{\Omega}_2 = - xy`                                                                  |
|             |                                  |             +----------------+------------------------------------------------------------------------------------------------+
|             |                                  |             | XZSECMOM       | :math:`\hat{\Omega}_3 = - xz`                                                                  |
|             |                                  |             +----------------+------------------------------------------------------------------------------------------------+
|             |                                  |             | YYSECMOM       | :math:`\hat{\Omega}_4 = - y^2`                                                                 |
|             |                                  |             +----------------+------------------------------------------------------------------------------------------------+
|             |                                  |             | YZSECMOM       | :math:`\hat{\Omega}_5 = - yz`                                                                  |
|             |                                  |             +----------------+------------------------------------------------------------------------------------------------+
|             |                                  |             | ZZSECMOM       | :math:`\hat{\Omega}_6 = - z^2`                                                                 |
+-------------+----------------------------------+-------------+----------------+------------------------------------------------------------------------------------------------+
| THETA       | Traceless quadrupole             |     S       | XXTHETA        | :math:`\hat{\Omega}_1 = - \frac{3}{2} (x^2-\frac{1}{3}r^2)`                                    |
|             |                                  |             +----------------+------------------------------------------------------------------------------------------------+
|             |                                  |             | XYTHETA        | :math:`\hat{\Omega}_2 = - \frac{3}{2} xy`                                                      |
|             |                                  |             +----------------+------------------------------------------------------------------------------------------------+
|             |                                  |             | XZTHETA        | :math:`\hat{\Omega}_3 = - \frac{3}{2} xz`                                                      |
|             |                                  |             +----------------+------------------------------------------------------------------------------------------------+
|             |                                  |             | YYTHETA        | :math:`\hat{\Omega}_4 = - \frac{3}{2} (y^2-\frac{1}{3}r^2)`                                    |
|             |                                  |             +----------------+------------------------------------------------------------------------------------------------+
|             |                                  |             | YZTHETA        | :math:`\hat{\Omega}_5 = - \frac{3}{2} yz`                                                      |
|             |                                  |             +----------------+------------------------------------------------------------------------------------------------+
|             |                                  |             | ZZTHETA        | :math:`\hat{\Omega}_6 = - \frac{3}{2} (z^2-\frac{1}{3}r^2)`                                    |
+-------------+----------------------------------+-------------+----------------+------------------------------------------------------------------------------------------------+
| SEFGMG      | Magnetic term of Flambaum-Ginges |     A       | XSEFGMG        |                                                                                                |
|             | self-energy potential            |             +----------------+------------------------------------------------------------------------------------------------+
|             |                                  |             | YSEFGMG        |                                                                                                |
|             |                                  |             +----------------+------------------------------------------------------------------------------------------------+
|             |                                  |             | ZSEFGMG        |                                                                                                |
+-------------+----------------------------------+-------------+----------------+------------------------------------------------------------------------------------------------+

where :math:`\boldsymbol{R}_{CM}` and :math:`\boldsymbol{R}_K` are the position vectors of the nuclear center of mass and the nucleus :math:`K` with respect to the origin of coordinates, respectively. In addition, the operators :math:`x`, :math:`y`, and :math:`z` are the Cartesian components of the electron position vector operator with respect to the origin of coordinates.

Note that, in the table above, :math:`g_e` is the electronic g-factor (i.e., :math:`g_e = 2`), 'NAM' are the first three letters of the name of atom :math:`K` as given in the .mol file, 'ab' denotes the two-digit index of the symmetry-adapted nucleus :math:`K`, and 'lmn', 'opq', and 'rst' denote three-digit indices identifying, respectively, the :math:`x`, :math:`y`, and :math:`z` Cartesian components of the given operator associated with nucleus :math:`K` in the symmetry-independent (non-symmetry-adapted) representation.

For the normalized Gaussian nuclear density distribution, the default exponents (:math:`\zeta_K`) are the values proposed by Visscher and Dyall in :cite:`Visscher1997b`.

The matrix representations of the operators shown in the previous Table are either symmetric (S) or antisymmetric (A). This information is provided in the third column of the Table.


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
UEHLING     Uehling potential
SEFGLF      Low-frequency term of Flambaum-Ginges self-energy potential
SEFGHF      High-frequency term of Flambaum-Ginges self-energy potential
SEPZ        Pyykko-Zhao self-energy model potential
=========== ====================================================================================================================


Examples of using various operators
-----------------------------------

We give here several concrete examples on how to construct operators for various properties.


Kinetic part of the Dirac Hamiltonian
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The kinetic part of the electronic Dirac Hamiltonian is given (in a.u.) by

.. math::

  \hat{H}^K & = c\vec{\boldsymbol{\alpha}}\cdot\hat{\vec{p}} = c\left(\boldsymbol{\alpha}_x\hat{p}_x+\boldsymbol{\alpha}_y\hat{p}_y+\boldsymbol{\alpha}_z\hat{p}_z\right)

            & = -c \left[ i\boldsymbol{\alpha}_x \frac{\partial}{\partial x} + i\boldsymbol{\alpha}_y \frac{\partial}{\partial y} + i\boldsymbol{\alpha}_z \frac{\partial}{\partial z} \right].

This operator may be specified by::

  .OPERATOR
   'Kin energy'
   ALPHADOT
   XDIPVEL
   YDIPVEL
   ZDIPVEL
   COMFACTOR
   -137.03599926085

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
   -1


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

The FC integrals are called using the label ``'FC NAMab'`` (see `One-electron scalar operator labels`_), and the total perturbed Hamiltonian :math:`\hat{H}=\hat{H}_0+\hat{H}'` (where, in a.u., :math:`\hat{H}'=10^{-9} \, \frac{4\pi}{3} \, g_e \, \delta(\boldsymbol{r}-\boldsymbol{R}_\text{Pb})`) can be requested using::

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

The sixteen :math:`4 \times 4` matrices listed in Section `Syntax for specifying one-electron operators`_, i.e., :math:`\boldsymbol{I}`, :math:`\boldsymbol{\beta}`, :math:`\boldsymbol{\gamma}^5`, :math:`\boldsymbol{\beta}\boldsymbol{\gamma}^5`, :math:`i\boldsymbol{\alpha}_x`, :math:`i\boldsymbol{\alpha}_y`, :math:`i\boldsymbol{\alpha}_z`, :math:`i\boldsymbol{\beta}\boldsymbol{\alpha}_x`, :math:`i\boldsymbol{\beta}\boldsymbol{\alpha}_y`, :math:`i\boldsymbol{\beta}\boldsymbol{\alpha}_z`, :math:`i\boldsymbol{\Sigma}_x`, :math:`i\boldsymbol{\Sigma}_y`, :math:`i\boldsymbol{\Sigma}_z`, :math:`i\boldsymbol{\beta}\boldsymbol{\Sigma}_x`, :math:`i\boldsymbol{\beta}\boldsymbol{\Sigma}_y`, and :math:`i\boldsymbol{\beta}\boldsymbol{\Sigma}_z`, are usually expressed in Molecular Physics using the standard representation where, for example,

.. math::

   \tilde{\boldsymbol{\beta}} =
   \begin{pmatrix}
     1 & 0 & 0  &  0 \\
     0 & 1 & 0  &  0 \\
     0 & 0 & -1 &  0 \\
     0 & 0 & 0  & -1
   \end{pmatrix}, \;
   \tilde{\boldsymbol{\gamma}}^{5} =
   \begin{pmatrix}
     0 & 0 & 1 & 0 \\
     0 & 0 & 0 & 1 \\
     1 & 0 & 0 & 0 \\
     0 & 1 & 0 & 0
   \end{pmatrix}, \;
   i\tilde{\boldsymbol{\alpha}}_z =
   \begin{pmatrix}
     0  & 0  & i  &  0 \\
     0  & 0  & 0  & -i \\
     i  & 0  & 0  &  0 \\
     0  & -i & 0  &  0
   \end{pmatrix}.

Here and in what follows, we use the notation :math:`\tilde{\boldsymbol{U}}` to refer to any of these sixteen matrices when they are expressed in the standard representation.

It can be seen that all these matrices can be written as the Kronecker product of two :math:`2 \times 2` matrices :math:`\boldsymbol{A}` and :math:`\boldsymbol{B}`, such that :math:`\tilde{\boldsymbol{U}} = \boldsymbol{A} \otimes \boldsymbol{B}`. In particular, :math:`\boldsymbol{A}` is one of the real matrices :math:`\boldsymbol{I}`, :math:`\boldsymbol{\sigma}_x`, :math:`i \boldsymbol{\sigma}_y`, or :math:`\boldsymbol{\sigma}_z`, whereas :math:`\boldsymbol{B}` is :math:`\boldsymbol{I}`, :math:`i \boldsymbol{\sigma}_x`, :math:`i \boldsymbol{\sigma}_y`, or :math:`i \boldsymbol{\sigma}_z`. Here, :math:`\boldsymbol{\sigma}_x`, :math:`\boldsymbol{\sigma}_y`, and :math:`\boldsymbol{\sigma}_z` are the Pauli matrices, and :math:`i` is the imaginary unit.

The individual matrix elements of :math:`\tilde{\boldsymbol{U}}` are given as

.. math::

    \tilde{U}_{pq} = A_{ij} \, B_{kl} \, ,

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

The same reordering applies to all the matrices :math:`\tilde{\boldsymbol{U}}`, which are then transformed accordingly. In the original ordering (i.e., in the standard representation), these matrices are given as

.. math::

    \tilde{\boldsymbol{U}} =
    \begin{pmatrix}
    \tilde{U}_{L\alpha,L\alpha} & \tilde{U}_{L\alpha,L\beta} & \tilde{U}_{L\alpha,S\alpha} & \tilde{U}_{L\alpha,S\beta} \\
    \tilde{U}_{L\beta,L\alpha}  & \tilde{U}_{L\beta,L\beta}  & \tilde{U}_{L\beta,S\alpha}  & \tilde{U}_{L\beta,S\beta}  \\
    \tilde{U}_{S\alpha,L\alpha} & \tilde{U}_{S\alpha,L\beta} & \tilde{U}_{S\alpha,S\alpha} & \tilde{U}_{S\alpha,S\beta} \\
    \tilde{U}_{S\beta,L\alpha}  & \tilde{U}_{S\beta,L\beta}  & \tilde{U}_{S\beta,S\alpha}  & \tilde{U}_{S\beta,S\beta}
    \end{pmatrix}
    = \boldsymbol{A} \otimes \boldsymbol{B}.

However, after reordering one clearly has the transformed matrices

.. math::

    \boldsymbol{U} =
    \begin{pmatrix}
    \tilde{U}_{L\alpha,L\alpha} & \tilde{U}_{L\alpha,S\alpha} & \tilde{U}_{L\alpha,L\beta} & \tilde{U}_{L\alpha,S\beta} \\
    \tilde{U}_{S\alpha,L\alpha} & \tilde{U}_{S\alpha,S\alpha} & \tilde{U}_{S\alpha,L\beta} & \tilde{U}_{S\alpha,S\beta} \\
    \tilde{U}_{L\beta,L\alpha}  & \tilde{U}_{L\beta,S\alpha}  & \tilde{U}_{L\beta,L\beta}  & \tilde{U}_{L\beta,S\beta}  \\
    \tilde{U}_{S\beta,L\alpha}  & \tilde{U}_{S\beta,S\alpha}  & \tilde{U}_{S\beta,L\beta}  & \tilde{U}_{S\beta,S\beta}
    \end{pmatrix},

which satisfy

.. math::

    \boldsymbol{U} = \boldsymbol{P} \, \tilde{\boldsymbol{U}} \, \boldsymbol{P}^{-1} =
    \boldsymbol{P} \, ( \boldsymbol{A} \otimes \boldsymbol{B} ) \, \boldsymbol{P}^{-1} =
    \boldsymbol{B} \otimes \boldsymbol{A},

where the :math:`4 \times 4` Kronecker permutation matrix :math:`\boldsymbol{P}` is a special permutation matrix that swaps the order of the Kronecker product of two :math:`2 \times 2` matrices, and is given by

.. math::

   \boldsymbol{P} =
   \begin{pmatrix}
   1 & 0 & 0 & 0 \\
   0 & 0 & 1 & 0 \\
   0 & 1 & 0 & 0 \\
   0 & 0 & 0 & 1
   \end{pmatrix}.

[Note that applying :math:`\boldsymbol{P}` on the left permutes the rows, while applying :math:`\boldsymbol{P}^{-1}` on the right permutes the columns, ensuring the correct transformation from :math:`\boldsymbol{A} \otimes \boldsymbol{B}` to :math:`\boldsymbol{B} \otimes \boldsymbol{A}`.]

In other words, :math:`U_{PQ} = B_{kl} \, A_{ij}`, with :math:`P=2(k-1)+i` and :math:`Q=2(l-1)+j`. This reordering guarantees that all the matrices :math:`\boldsymbol{U}` are time-reversal symmetric.


Additionally, following the arguments given in Ref. :cite:`Saue1999`, it is possible to make the following mapping involving the quaternion units:

.. math::

    &\boldsymbol{I}          &\rightarrow & \qquad  1

    &i \boldsymbol{\sigma}_z &\rightarrow & \qquad  \check{\imath}

    &i \boldsymbol{\sigma}_y &\rightarrow & \qquad  \check{\jmath}

    &i \boldsymbol{\sigma}_x &\rightarrow & \qquad  \check{k}

or, in a general way, the :math:`2 \times 2` matrices :math:`\boldsymbol{B}` can be mapped as :math:`\boldsymbol{B} \rightarrow b`, where :math:`b` is a quaternion unit (i.e., :math:`1`, :math:`\check{\imath}`, :math:`\check{\jmath}`, or :math:`\check{k}`) arising from the relationships given above.

Thus, any linear combination of the 16 reordered matrices :math:`\boldsymbol{U}` can be mapped to a linear combination of :math:`2 \times 2` matrices :math:`\boldsymbol{u}`, each of which contains a quaternion unit.


To make this discussion clearer, in the following we discuss two examples. In first place, we transform the :math:`4 \times 4` matrix :math:`\tilde{\boldsymbol{\beta}}`. We start by writting this matrix as the following Kronecker product of two :math:`2 \times 2` matrices: :math:`\tilde{\boldsymbol{\beta}}=\boldsymbol{\sigma}_z \otimes \boldsymbol{I}`. In the reordered representation, this matrix will be expressed as :math:`\boldsymbol{\beta}=\boldsymbol{I} \otimes \boldsymbol{\sigma}_z`. Finally, the corresponding :math:`\boldsymbol{B}` matrix (in this case, :math:`\boldsymbol{I}`) is transformed to a quaternion unit (in the current example, :math:`1`). In this way, the :math:`4 \times 4` reordered matrix :math:`\boldsymbol{\beta}` can be written in quaternion form as the :math:`2 \times 2` matrix :math:`\boldsymbol{\sigma}_z`. In other words:

.. math::

   \tilde{\boldsymbol{\beta}} =
   \begin{pmatrix}
     1 & 0 & 0  &  0 \\
     0 & 1 & 0  &  0 \\
     0 & 0 & -1 &  0 \\
     0 & 0 & 0  & -1
   \end{pmatrix}
   \;\rightarrow\;
   \boldsymbol{\beta} =
   \boldsymbol{P}\,\tilde{\boldsymbol{\beta}}\,\boldsymbol{P}^{-1} =
   \begin{pmatrix}
     1 & 0  & 0 &  0 \\
     0 & -1 & 0 &  0 \\
     0 & 0  & 1 &  0 \\
     0 & 0  & 0 & -1
   \end{pmatrix}
   \;\rightarrow\;
   \boldsymbol{\beta}^\text{quat} =
   \boldsymbol{\sigma}_z.


As a second example, we take :math:`i\tilde{\boldsymbol{\alpha}}_z`. This :math:`4 \times 4` matrix, in its original ordering (according to the standard representation), can be written as the Kronecker product :math:`i\tilde{\boldsymbol{\alpha}}_z=\boldsymbol{\sigma}_x \otimes i \boldsymbol{\sigma}_z`. This means that in the reordered representation it will be given as :math:`i\boldsymbol{\alpha}_z = i \boldsymbol{\sigma}_z \otimes \boldsymbol{\sigma}_x`. In this particular case, the corresponding matrix :math:`\boldsymbol{B}` is :math:`i \boldsymbol{\sigma}_z`, so that the quaternion unit :math:`b` is :math:`\check{\imath}`. Therefore, the quaternion form of this reordered matrix will be :math:`i\boldsymbol{\alpha}_z^\text{quat} = \check{\imath} \boldsymbol{\sigma}_x`. Equivalently:

.. math::

   i\tilde{\boldsymbol{\alpha}}_z =
   \begin{pmatrix}
     0 & 0  & i  &  0  \\
     0 & 0  & 0  &  -i \\
     i & 0  & 0  &  0  \\
     0 & -i & 0  &  0
   \end{pmatrix}
   \;\rightarrow\;
   i\boldsymbol{\alpha}_z =
   \boldsymbol{P}\,i\tilde{\boldsymbol{\alpha}}_z\,\boldsymbol{P}^{-1} =
   \begin{pmatrix}
     0 & i  & 0  &  0  \\
     i & 0  & 0  &  0  \\
     0 & 0  & 0  &  -i \\
     0 & 0  & -i &  0
   \end{pmatrix}
   \;\rightarrow\;
   i\boldsymbol{\alpha}_z^\text{quat} =
   \check{\imath} \boldsymbol{\sigma}_x.


The relationships between matrices :math:`\tilde{\boldsymbol{U}}` (i.e., the sixteen :math:`4 \times 4` matrices in original ordering, which corresponds to the standard representation) and the quaternion form of their reordered equivalent :math:`2 \times 2` matrices :math:`\boldsymbol{u}` is given in the following Table. They follow the transformation rule :math:`\tilde{\boldsymbol{U}} = \boldsymbol{A} \otimes \boldsymbol{B}` :math:`\rightarrow` :math:`\boldsymbol{U} = \boldsymbol{P} \, \tilde{\boldsymbol{U}} \, \boldsymbol{P}^{-1} = \boldsymbol{P} \, (\boldsymbol{A} \otimes \boldsymbol{B}) \, \boldsymbol{P}^{-1} = \boldsymbol{B} \otimes \boldsymbol{A}` :math:`\rightarrow` :math:`\boldsymbol{u} = b \boldsymbol{A}`.


+------------------------------------------------------------------+------------------------------------------------------------------+---------------------------------------------------+
| :math:`\tilde{\boldsymbol{U}}`                                   | :math:`\boldsymbol{A} \otimes \boldsymbol{B}`                    | :math:`\boldsymbol{u}`                            |
+==================================================================+==================================================================+===================================================+
| :math:`\tilde{\boldsymbol{I}}`                                   | :math:`\boldsymbol{I} \otimes \boldsymbol{I}`                    | :math:`\boldsymbol{I}`                            |
+------------------------------------------------------------------+------------------------------------------------------------------+---------------------------------------------------+
| :math:`i\tilde{\boldsymbol{\Sigma}}_z`                           | :math:`\boldsymbol{I} \otimes i \boldsymbol{\sigma}_z`           | :math:`\check{\imath}\, \boldsymbol{I}`           |
+------------------------------------------------------------------+------------------------------------------------------------------+---------------------------------------------------+
| :math:`i\tilde{\boldsymbol{\Sigma}}_y`                           | :math:`\boldsymbol{I} \otimes i \boldsymbol{\sigma}_y`           | :math:`\check{\jmath}\, \boldsymbol{I}`           |
+------------------------------------------------------------------+------------------------------------------------------------------+---------------------------------------------------+
| :math:`i\tilde{\boldsymbol{\Sigma}}_x`                           | :math:`\boldsymbol{I} \otimes i \boldsymbol{\sigma}_x`           | :math:`\check{k}\, \boldsymbol{I}`                |
+------------------------------------------------------------------+------------------------------------------------------------------+---------------------------------------------------+
| :math:`\tilde{\boldsymbol{\beta}}`                               | :math:`\boldsymbol{\sigma}_z \otimes \boldsymbol{I}`             | :math:`\boldsymbol{\sigma}_z`                     |
+------------------------------------------------------------------+------------------------------------------------------------------+---------------------------------------------------+
| :math:`i\tilde{\boldsymbol{\beta}}\tilde{\boldsymbol{\Sigma}}_z` | :math:`\boldsymbol{\sigma}_z \otimes i \boldsymbol{\sigma}_z`    | :math:`\check{\imath}\, \boldsymbol{\sigma}_z`    |
+------------------------------------------------------------------+------------------------------------------------------------------+---------------------------------------------------+
| :math:`i\tilde{\boldsymbol{\beta}}\tilde{\boldsymbol{\Sigma}}_y` | :math:`\boldsymbol{\sigma}_z \otimes i \boldsymbol{\sigma}_y`    | :math:`\check{\jmath}\, \boldsymbol{\sigma}_z`    |
+------------------------------------------------------------------+------------------------------------------------------------------+---------------------------------------------------+
| :math:`i\tilde{\boldsymbol{\beta}}\tilde{\boldsymbol{\Sigma}}_x` | :math:`\boldsymbol{\sigma}_z \otimes i \boldsymbol{\sigma}_x`    | :math:`\check{k}\, \boldsymbol{\sigma}_z`         |
+------------------------------------------------------------------+------------------------------------------------------------------+---------------------------------------------------+
| :math:`\tilde{\boldsymbol{\beta}}\tilde{\boldsymbol{\gamma}}^5`  | :math:`(i\boldsymbol{\sigma}_y) \otimes \boldsymbol{I}`          | :math:`(i\boldsymbol{\sigma}_y)`                  |
+------------------------------------------------------------------+------------------------------------------------------------------+---------------------------------------------------+
| :math:`i\tilde{\boldsymbol{\beta}}\tilde{\boldsymbol{\alpha}}_z` | :math:`(i\boldsymbol{\sigma}_y) \otimes i \boldsymbol{\sigma}_z` | :math:`\check{\imath}\, (i\boldsymbol{\sigma}_y)` |
+------------------------------------------------------------------+------------------------------------------------------------------+---------------------------------------------------+
| :math:`i\tilde{\boldsymbol{\beta}}\tilde{\boldsymbol{\alpha}}_y` | :math:`(i\boldsymbol{\sigma}_y) \otimes i \boldsymbol{\sigma}_y` | :math:`\check{\jmath}\, (i\boldsymbol{\sigma}_y)` |
+------------------------------------------------------------------+------------------------------------------------------------------+---------------------------------------------------+
| :math:`i\tilde{\boldsymbol{\beta}}\tilde{\boldsymbol{\alpha}}_x` | :math:`(i\boldsymbol{\sigma}_y) \otimes i \boldsymbol{\sigma}_x` | :math:`\check{k}\, (i\boldsymbol{\sigma}_y)`      |
+------------------------------------------------------------------+------------------------------------------------------------------+---------------------------------------------------+
| :math:`\tilde{\boldsymbol{\gamma}}^5`                            | :math:`\boldsymbol{\sigma}_x \otimes \boldsymbol{I}`             | :math:`\boldsymbol{\sigma}_x`                     |
+------------------------------------------------------------------+------------------------------------------------------------------+---------------------------------------------------+
| :math:`i\tilde{\boldsymbol{\alpha}}_z`                           | :math:`\boldsymbol{\sigma}_x \otimes i \boldsymbol{\sigma}_z`    | :math:`\check{\imath}\, \boldsymbol{\sigma}_x`    |
+------------------------------------------------------------------+------------------------------------------------------------------+---------------------------------------------------+
| :math:`i\tilde{\boldsymbol{\alpha}}_y`                           | :math:`\boldsymbol{\sigma}_x \otimes i \boldsymbol{\sigma}_y`    | :math:`\check{\jmath}\, \boldsymbol{\sigma}_x`    |
+------------------------------------------------------------------+------------------------------------------------------------------+---------------------------------------------------+
| :math:`i\tilde{\boldsymbol{\alpha}}_x`                           | :math:`\boldsymbol{\sigma}_x \otimes i \boldsymbol{\sigma}_x`    | :math:`\check{k}\, \boldsymbol{\sigma}_x`         |
+------------------------------------------------------------------+------------------------------------------------------------------+---------------------------------------------------+


It is important to note that although the time-reversal antisymmetric operators :math:`\boldsymbol{\alpha}_k`, :math:`\boldsymbol{\Sigma}_k`, :math:`\boldsymbol{\beta\alpha}_k`, and :math:`\boldsymbol{\beta\Sigma}_k` (with :math:`k=x,y`, and :math:`z`) become time-reversal symmetric once they are multiplied by the imaginary unit :math:`i` and written in quaternion form as described above, DIRAC saves the information that these operators were originally T-reversal antisymmetric (before adding the imaginary phase).

For further details on the quaternion symmetry scheme used in DIRAC, see Ref. :cite:`Saue1999`.

`Note:` For complete consistency, the names of the operators related to :math:`i\boldsymbol{\alpha}_k`, :math:`i\boldsymbol{\Sigma}_k`, :math:`i\boldsymbol{\beta\alpha}_k`, and :math:`i\boldsymbol{\beta\Sigma}_k` (:math:`k=x,y,z`) could have been ``XiALPHA``, ``YiALPHA``, ``ZiALPHA``, ``XiAVECTOR``, ``YiAVECTOR``, ``ZiAVECTOR``, ``iALPHADOT``, ``XiBETAALP``, ``YiBETAALP``, ``ZiBETAALP``, ``XiSIGMA``, ``YiSIGMA``, ``ZiSIGMA``, ``XiBETASIG``, ``YiBETASIG``, ``ZiBETASIG``, ``iSIGMADOT``, and ``iBEALDOT`` (see `Operator type keywords`_). However, due to historical reasons, this convention is not used in DIRAC.

`Note for advanced users:` The operator iBETAGAM is implemented in the code just as :math:`\beta\boldsymbol{\gamma}^5`, but its time-reversal symmetry without imaginary phase is set as odd.
