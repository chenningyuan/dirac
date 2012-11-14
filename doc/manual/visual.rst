

.. _visual:
.. index:: **VISUAL
==========
\*\*VISUAL
==========


.. index:: .LIST
.LIST
-----

Calculate various densities for
a list of points. Example (3 points; coordinates in bohr)::

  .LIST
   3
   1.0 0.0 0.0
   0.0 1.0 0.0
   0.0 0.0 1.0


.. index:: .LINE
.LINE
-----

Calculate various densities along a line.
Example (line connecting two points; 200 steps; coordinates in bohr)::

  .LINE
   0.0 0.0 0.0
   0.0 0.0 5.0
   200


.. index:: .2D
.2D
---

Calculate various densities in a plane.
The plane is specified using 3 points that have to form a right angle.
Example (coordinates in bohr)::

  .2D
   0.0  0.0  0.0     !origin
   0.0  0.0 10.0     !"right"
   200               !nr of points origin-"right"
   0.0 10.0  0.0     !"top"
   200               !nr of points origin-"top"


.. index:: .2DINT
.2DINT
------

Integrate various densities in a plane.
The plane is specified using 3 points that have to form a right angle.
Example (coordinates in bohr)::

  .2DINT
   0.0  0.0  0.0     !origin
   0.0  0.0 10.0     !"right"
   200               !nr of points origin-"right"
   0.0 10.0  0.0     !"top"
   200               !nr of points origin-"top"


.. index:: .3D
.3D
---

Calculate various densities in 3D and write to cube file format.
Example (coordinates in bohr)::

  .3D
   40 40 40          ! 40 x 40 x 40 points


.. index:: .3DFAST
.3DFAST
-------

Fast evaluation of the molecular electrostatic potential.
Example (coordinates in bohr)::

  .3DFAST
   40 40 40          ! 40 x 40 x 40 points


.. index:: .3DADD
.3DADD
------

Add space around the cube file
Default (coordinates in bohr)::

  .3D_ADD
   4.0


.. index:: .3D_INT
.3D_INT
-------

Integrate densities in 3D.


.. index:: .SCALE
.SCALE
------

Scale densities by a factor.
Default::

  .SCALE
   1.0


.. index:: .DENSITY
.DENSITY
--------

Compute density.
Example (unperturbed density)::

  .DENSITY
   DFCOEF

Another example (perturbed density, first response vector)::

  .DENSITY
   PAMXVC 1


.. index:: .ELF
.ELF
----

Compute the electron localization function. Example::

  .ELF
   DFCOEF


.. index:: .GAMMA5
.GAMMA5
-------

Compute the electron chirality density. Example::

  .GAMMA5
   DFCOEF


.. index:: .JDIA
.JDIA
-----

Compute the nonrelativistic diamagnetic current density. Example::

  .JDIA
   DFCOEF


.. index:: .J
.J
--

Compute the current density. Example (use first response vector)::

  .J
   PAMXVC 1


.. index:: .DIVJ
.DIVJ
-----

Compute the divergence of the current density. Example (use first response vector)::

  .DIVJ
   PAMXVC 1


.. index:: .ROTJ
.ROTJ
-----

Compute the curl of the current density. Example (use first response vector)::

  .ROTJ
   PAMXVC 1


.. index:: .ESP
.ESP
----

Compute the electrostatic potential. Example::

  .ESP
   DFCOEF


.. index:: .ESPE
.ESPE
-----

Compute the electronic part of the electrostatic potential.


.. index:: .ESPE
.ESPN
-----

Compute the nuclear part of the electrostatic potential.


.. index:: .ESPRHO
.ESPRHO
-------

Compute the electrostatic potential times density.


.. index:: .ESPERHO
.ESEPRHO
--------

Compute the electronic part of the electrostatic potential times density.


.. index:: .ESPNRHO
.ESNPRHO
--------

Compute the nuclear part of the electrostatic potential times density.


.. index:: .GAUGE
.GAUGE
------

Specify gauge origin. Example::

  .GAUGE
   0.0 0.0 0.0


.. index:: .RADIAL
.RADIAL
-------

Compute radial densities.


.. index:: .SMALLAO
.SMALLAO
--------

Force evaluation of small component basis functions.


.. index:: .OCCUPATION
.OCCUPATION
-----------

Specify occupation of orbitals.
Example (neon atom)::

  .OCCUPATION
   2
   1 1-2 1.0
   2 1-3 1.0

Another example (water)::

  .OCCUPATION
   1
   1 1-5 1.0

Another example (nitrogen atom)::

  .OCCUPATION
   2
   1 1-2 1.0
   2 1-3 0.5

Using dominating orbital rotation amplitudes you can use this to plot
transition densities.


.. index:: .LONDON
.LONDON
-------

Activate LAO contribution.


.. index:: .NONE
.NONE
-----

Select "none" connection when
when plotting LAO perturbed densities.


.. index:: .NODIRECT
.NODIRECT
---------

Skip direct LAO contribution
when plotting perturbed densities.


.. index:: .NOREORTHO
.NOREORTHO
----------

Skip LAO reorthonormalization contribution
when plotting perturbed densities.


.. index:: .NOKAPPA
.NOKAPPA
--------

Skip orbital relaxation contribution
when plotting perturbed densities.
