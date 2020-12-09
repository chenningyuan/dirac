:orphan:
 

.. warning::

   Only the calculation of the density is tested for open shell configurations
   (and relies on a correct .OCCUPATION). All other properties are only
   tested for closed shell systems and should not be trusted for open shell systems
   without a thorough testing.


starstar(VISUAL)

Sampling
========

keyword(LIST)

Calculate various densities in few points. Scalar and vector densities are written to the standard output file. Example (3 points; coordinates in
bohr)::

  .LIST
   3
   1.0 0.0 0.0
   0.0 1.0 0.0
   0.0 0.0 1.0


keyword(LINE)

Calculate various densities along a line.
Example (line connecting two points; 200 steps; coordinates in bohr)::

  .LINE
   0.0 0.0 0.0
   0.0 0.0 5.0
   200

Scalar and vector densities are written to files ``plot.line.scalar`` and ``plot.line.vector``, respectively, 
and should be saved after calculation, e.g. ::

   pam --get=plot.line.scalar ...

The first three columns of the output files gives the coordinates (x, y, z) of the point. 
It is then followed by one/three columns giving the value of the scalar/vector density in that point.

keyword(RADIAL)

Compute radial distributions

.. math::

   f(r) = \int_{0}^{2\pi}\int_{0}^{\pi}f(\mathbf{r})r^2\sin\theta d\theta d\phi

by performing Lebedev angular integration over a specified number of even-spaced radial shells out to some specified distance
from a specified initial point. Example (coordinates and distance in bohr)::

   .RADIAL
   0.0 0.0 0.0
   10.0
   200

The first line after the keyword specifies the initial point, here chosen to be the origin.
The second and third line is the distance and step size, respectively.
Scalar and vector densities are written to files ``plot.radial.scalar`` and ``plot.radial.vector``, respectively, 
and should be saved after calculation, e.g. ::

   pam --get=plot.radial.scalar ...

keyword(2D)

Calculate various densities in a plane.
The plane is specified using 3 points that have to form a right angle.
Example (coordinates in bohr)::

  .2D
   0.0  0.0  0.0     !origin
   0.0  0.0 10.0     !"right"
   200               !nr of points origin-"right"
   0.0 10.0  0.0     !"top"
   200               !nr of points origin-"top"

Scalar and vector densities are written to files ``plot.2d.scalar`` and ``plot.2d.vector``, respectively, and should be saved after calculation, e.g. ::

   pam --get=plot.2d.scalar ...

keyword(2D_INT)

Integrate various densities in a plane
using Gauss-Lobatto quadrature.
The plane is specified using 3 points that have to form a right angle.
Example (coordinates in bohr)::

  .2D_INT
   0.0  0.0  0.0     !origin
   0.0  0.0 10.0     !"right"
   10                !nr of tiles to the "right"
   0.0 10.0  0.0     !"top"
   10                !nr of tiles to the "right"
   5                 !order of the Legendre polynomial for each tile


keyword(3D)

Calculate various densities in 3D and write to cube file format.
Example (coordinates in bohr)::

  .3D
   40 40 40          ! 40 x 40 x 40 points

Scalar and vector densities are written to files ``plot.3d.scalar`` and ``plot.3d.vector``, respectively, and should be saved after calculation, e.g. ::

   pam --get=plot.3d.scalar ...
   
Scalar densities are also written to a Gaussian cube file ``plot.3d.cube``.

keyword(3D_IMP)

Calculate various densities in 3D on an imported grid and write to text file format (imported grid does not have to be regular).
Example:

  .3D_IMP
   numerical_grid    ! a grid file on which various densities will be calculated; the file should have 3 columns with x, y, z coordinates of grid points 
   grid_units        ! grid units (choose from: 'angstrom' or 'au')


keyword(3DFAST)

Fast evaluation of the molecular electrostatic potential.
Example (coordinates in bohr)::

  .3DFAST
   40 40 40          ! 40 x 40 x 40 points


keyword(3D_ADD)

Add space around the cube file.
Default (coordinates in bohr)::

  .3D_ADD
   4.0


keyword(3D_IMP)

Calculate various densities in 3D on an imported grid (does not have to be regular)
Example::

  .3D_IMP
   grid_file        ! a file with x,y,z-coordinates of grid points


keyword(3D_INT)

Integrate densities in 3D.

Modification of densities
=========================

keyword(CARPOW)

Scale densities by Cartesian product :math:`x^iy^jz^k`. The keyword is followed by three integers specifying the exponents :math:`(i,j,k)`. Example::

  .DENSITY
  .CARPOW
  1 0 0

is equivalent to the specification::

   .EDIPX

keyword(SCALE)

Scale densities by a factor.
Default::

  .SCALE
   1.0

keyword(DSCALE)

Scale densities *down* by a factor.
Default::

  .DSCALE
   1.0

Densities
=========   
   
keyword(DENSITY)

Compute number density :math:`n(\mathbf{r})` .
Example (unperturbed density)::

  .DENSITY
   DFCOEF

Another example (perturbed density, first response vector)::

  .DENSITY
   PAMXVC 1

   
keyword(DENGRA)

Compute number density gradient :math:`\grad n(\mathbf{r})` in addition to the number density :math:`n(\mathbf{r})`.
Example (unperturbed density gradient)::

  .DENGRA
   DFCOEF


keyword(ELF)

Compute the electron localization function. Example::

  .ELF
   DFCOEF


keyword(GAMMA5)

Compute the electron chirality density. Example::

  .GAMMA5
   DFCOEF


keyword(J)

Compute the current density :math:`\mathbf{j}(\mathbf{r})=-e\psi_{i}^{\ast}c\boldsymbol{\alpha}\psi_{i}`. Example (use first response vector)::

  .J
   PAMXVC 1

keyword(JDIA)

Compute the nonrelativistic diamagnetic current density. Example::

  .JDIA
   DFCOEF


keyword(JX)

Compute the x-component :math:`j_{x}(\mathbf{r})=-e\psi_{i}^{\ast}c\alpha_{x}\psi_{i}` of the current density. Example (use first response vector)::

  .JX
   PAMXVC 1

keyword(JY)

Compute the y-component :math:`j_{y}(\mathbf{r})=-e\psi_{i}^{\ast}c\alpha_{y}\psi_{i}` of the current density. Example (use first response vector)::

  .JY
   PAMXVC 1

keyword(JZ)

Compute the z-component :math:`j_{z}(\mathbf{r})=-e\psi_{i}^{\ast}c\alpha_{z}\psi_{i}` of the current density. Example (use first response vector)::

  .JZ
   PAMXVC 1

keyword(DIVJ)

Compute the divergence of the current density. Example (use first response vector)::

  .DIVJ
   PAMXVC 1


keyword(ROTJ)

Compute the curl of the current density. Example (use first response vector)::

  .ROTJ
   PAMXVC 1

keyword(BDIPX)

Compute the x-component :math:`m^{[1]}_{x}(\mathbf{r})=-\frac{1}{2}(\mathbf{r}\times\mathbf{j})_{x}` of the magnetic dipole operator. Example (use first response vector)::

  .BDIPX
   PAMXVC 1

keyword(BDIPY)

Compute the y-component :math:`m^{[1]}_{y}(\mathbf{r})=-\frac{1}{2}(\mathbf{r}\times\mathbf{j})_{y}` of the magnetic dipole operator. Example (use first response vector)::

  .BDIPY
   PAMXVC 1

keyword(BDIPZ)

Compute the z-component :math:`m^{[1]}_{z}(\mathbf{r})=-\frac{1}{2}(\mathbf{r}\times\mathbf{j})_{z}` of the magnetic dipole operator. Example (use first response vector)::

  .BDIPZ
   PAMXVC 1

keyword(EDIPX)

Compute the x-component :math:`Q^{[1]}_{x}(\mathbf{r})=xn(\mathbf{r})` of the electric dipole.

keyword(EDIPY)

Compute the y-component :math:`Q^{[1]}_{y}(\mathbf{r})=yn(\mathbf{r})` of the electric dipole.

keyword(EDIPZ)

Compute the z-component :math:`Q^{[1]}_{z}(\mathbf{r})=zn(\mathbf{r})` of the electric dipole.

keyword(ESP)

Compute the electrostatic potential. Example::

  .ESP
   DFCOEF


keyword(ESPE)

Compute the electronic part of the electrostatic potential.


keyword(ESPN)

Compute the nuclear part of the electrostatic potential.


keyword(ESPRHO)

Compute the electrostatic potential times density.


keyword(ESPERHO)

Compute the electronic part of the electrostatic potential times density.


keyword(ESPNRHO)

Compute the nuclear part of the electrostatic potential times density.


keyword(NDIPX)

Compute the NMR shielding density, with the "X"-component of the nuclear magnetic dipole moment and
the selected component of the magnetically-induced current density (by the chosen record on PAMXVC file) as perturbing operators.

keyword(NDIPY)

Compute the NMR shielding density, with the "Y"-component of the nuclear magnetic dipole moment and
the selected component of the magnetically-induced current density (by the chosen record on PAMXVC file) as perturbing operators.

keyword(NDIPZ)

Compute the NMR shielding density, with the "Z"-component of the nuclear magnetic dipole moment and
the selected component of the magnetically-induced current density (by the chosen record on PAMXVC file) as perturbing operators.

keyword(NICS )

Compute the NMR shielding density in a selected point in space. Is used to calculate NICS.
Example::

  .NICS
   1.2 -1.0 2.0

will calculate the NMR shielding in point (1.2, -1.0, 2.0). This keyword can be used only with one of: NDIPX, NDIPY, NDIPZ keywords.


keyword(READJB)

Use the grid and the magnetically-induced current density (jB) from a file to calculate the jB-dependent densities,
e.g. the NMR shielding density or the magnetizability density. 
Example::

  .READJB
   file_name        ! a file with x,y,z-coordinates of grid points and jB vector field


keyword(SMALLAO)

Force evaluation of small component basis functions.


keyword(OCCUPATION)

Specify occupation of orbitals.
Example (neon atom)::

  .OCCUPATION
   2
   1 1-2 1.0
   2 1-3 1.0

The first line after the keyword gives the number of subsequent lines to read.
In each line, the first number is the fermion ircop. In molecules with inversion symmetry 
there are two fermion ircops: gerade (1) and ungerade (2). Otherwise there is a single fermion ircop (1).
The specification of the fermion ircop is followed by the range of selected orbitals and their occupation.
If a single orbital is specified a single number is given instead of the range.

Another example (water)::

  .OCCUPATION
   1
   1 1-5 1.0

Another example (nitrogen atom)::

  .OCCUPATION
   2
   1 1-2 1.0
   2 1-3 0.5


keyword(LONDON)

Activate LAO contribution.
This keyword is followed by a letter "X", "Y" or "Z" indicating the component of an external perturbing magnetic field.
For example::

  .LONDON
   X


keyword(NONE)

Select "none" connection when
when plotting LAO perturbed densities.


keyword(NODIRECT)

Skip direct LAO contribution
when plotting perturbed densities.


keyword(NOREORTHO)

Skip LAO reorthonormalization contribution
when plotting perturbed densities.


keyword(NOKAPPA)

Skip orbital relaxation contribution
when plotting perturbed densities.
