

.. warning::

   Only the calculation of the density is tested for open shell configurations
   (and relies on a correct .OCCUPATION). All other properties are only
   tested for closed shell systems and should not be trusted for open shell systems
   without a thorough testing.


starstar(VISUAL)

keyword(LIST)

Calculate various densities in few points. Example (3 points; coordinates in
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


keyword(3DFAST)

Fast evaluation of the molecular electrostatic potential.
Example (coordinates in bohr)::

  .3DFAST
   40 40 40          ! 40 x 40 x 40 points


keyword(3DADD)

Add space around the cube file.
Default (coordinates in bohr)::

  .3D_ADD
   4.0


keyword(3D_INT)

Integrate densities in 3D.


keyword(SCALE)

Scale densities by a factor.
Default::

  .SCALE
   1.0


keyword(DENSITY)

Compute density.
Example (unperturbed density)::

  .DENSITY
   DFCOEF

Another example (perturbed density, first response vector)::

  .DENSITY
   PAMXVC 1


keyword(ELF)

Compute the electron localization function. Example::

  .ELF
   DFCOEF


keyword(GAMMA5)

Compute the electron chirality density. Example::

  .GAMMA5
   DFCOEF


keyword(JDIA)

Compute the nonrelativistic diamagnetic current density. Example::

  .JDIA
   DFCOEF


keyword(J)

Compute the current density. Example (use first response vector)::

  .J
   PAMXVC 1


keyword(DIVJ)

Compute the divergence of the current density. Example (use first response vector)::

  .DIVJ
   PAMXVC 1


keyword(ROTJ)

Compute the curl of the current density. Example (use first response vector)::

  .ROTJ
   PAMXVC 1


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


keyword(GAUGE)

Specify gauge origin. Example::

  .GAUGE
   0.0 0.0 0.0


keyword(RADIAL)

Compute radial densities.


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
