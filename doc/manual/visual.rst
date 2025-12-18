:orphan:
 
starstar(VISUAL)

This section describes the functionalities of the `VISUAL` module, which calculates various molecular densities from the 1-electron density matrices and atomic orbitals. These densities (here called `grid functions`) are calculated on numerical grids and can be integrated or exported as scalar, vector, or tensor fields (e.g., for visualization or further processing). The VISUAL module is not activated by default.

The main setup for these calculations should be provided by the `.GRIDS` and `.GRIDFUNCTIONS` blocks, as described below.

GRIDS
=====

keyword(GRIDS)

Grid(s) specification block. Next line should contain the number of grids to be used. Then, the following lines should provide information about each grid (as a `key=value` list).

For example, to import the grid from the `grid.h5` file, use::

  .GRIDS
  1
  id=1 input=import file_inp=grid.h5
  
Available `key=value` options
-----------------------------

* `id=int`

    The ID of a grid. Required, used as a reference for :ref:`VISUAL_.GRIDFUNCTIONS` specification.

* `npoints=[Nx,Ny,nz]`

   Specify the number of grid points in `x`, `y`, and `z` directions. Either `npoints` or `spacing` is required for grid generation.

* `spacing=[Sx,Sy,Sz]` or `spacing=S`

    Specify the grid resolution, i.e., the distance between two consecutive grid points. Use `spacing=[Sx,Sy,Sz]` to specify the distance between grid points in `x`, `y`, and `z` directions, or `spacing=S` to specify the same distance between grid points for all directions. Either `npoints` or `spacing` is required for grid generation.

* `input=option`

    Specify grid input. Available options:

        * `input=create` : create grid; the default is a 3D rectilinear grid; requires the specification of `npoints` or `spacing`
        * `input=import` : import grid from an external file
        * `input=dft` : use the DFT numerical grid; requires `numerical_grid` file
        * `input=list` : use the list of grid points; should be provided in the input (see :ref:`VISUAL.LIST`)
        * `input=line` : create grid points on line; should be provided in the input (see :ref:`VISUAL.LINE` keyword)
        * `input=2d` : create grid points on a 2D plane; should be provided in the input (see :ref:`VISUAL.2D` keyword)
        * `input=2d_int` : create grid points on a 2D plane for integration; should be provided in the input (see :ref:`VISUAL.2D_INT` keyword)
        * `input=3d` : create grid points on a 3D surface; should be provided in the input (see :ref:`VISUAL.3D` keyword)
        * `input=3d_int` : create grid points on a 3D surface for integration; should be provided in the input (see :ref:`VISUAL.3D_INT` keyword)
        * `input=radial` : create radial grid; should be provided in the input (see :ref:`VISUAL.RADIAL` keyword)

* `margin=f`

    Add space around the 3D rectilinear grid. Applies only to generated 3D rectilinear grids. Default is `f=0.0` (value in a.u.), which creates grid point spanning the space between the minimum and maximum coordinates of all atoms, i.e., :math:`x_{grid}=[min(x_{a}),max(x_{a})]_{a\in atoms}`, :math:`y_{grid}=[min(y_{a}),max(y_{a})]_{a\in atoms}`, :math:`z_{grid}=[min(z_{a}),max(z_{a})]_{a\in atoms}`. Therefore, for typical applications, choose some :math:`f>0` (e.g., `margin=4.0` is typically OK for the electron density visualization).

* `export=option`

    If the grid should be exported to a file, this allows to define the file format for export. Supported formats:

        * `export=hdf5` : export to hdf5 file
        * `export=txt` : export to text file, no header, space as a separator (for backward compatibility with old VISUAL code)
        * `export=csv` : export to text file (csv): header line, coma as a separator

* `file_inp=string`

    The name of a file, from which the grid should be read.

* `file_out=string`

    The name of a file, to which the grid should be exported

Notes:

* *current limitations:*
    * make sure there is one space between each `key=value` entry in the grid specification line.
    * the maximum line length in the input file is 200 characters
* for the specific examples, see the corresponding test in `DIRAC` test suite: `test/visual_grid_options`

Additional keywords
-------------------
The selected grid inputs (see `input=option`) require additional specification. This should be provided as a separate keyword block **following** the `.GRID` keyword. Most of these keywords are the same as used in the old VISUAL code.

keyword(LIST)

Calculate various densities in few points. Scalar and vector densities are written to the standard output file. Example (3 points; coordinates in bohr)::

  .LIST
   3
   1.0 0.0 0.0
   0.0 1.0 0.0
   0.0 0.0 1.0
   
In the `DIRAC` input file, this should be provided as::

  .GRIDS
  id=1 input=list
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

keyword(RADIAL)

Compute radial distributions

.. math::

   f(r) = \int_{0}^{2\pi}\int_{0}^{\pi}f(\mathbf{r})r^2\sin\theta d\theta d\phi

by performing Lebedev angular integration over a specified number of even-spaced radial shells out to some specified distance from a specified initial point. Example (coordinates and distance in bohr)::

   .RADIAL
   0.0 0.0 0.0
   10.0
   200

The first line after the keyword specifies the initial point, here chosen to be the origin.
The second and third line is the distance and step size, respectively.

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

Integrate various densities in a plane using Gauss-Lobatto quadrature.
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

Note: this is the same as using `npoints=[40,40,40]` in the grid specification line.


Grid functions
==============

keyword(GRIDFUNCTIONS) 

Grid function(s) specification block.
Next line should contain the number of grid functions to be calculated. Then, the following lines should provide information about each grid function (as a `key=value` list). 

For example, to calculate the values of the electron density (*ed*), its laplacian (*ed_laplacian*), and the reduced density gradient (*rdg*) in the grid points of the grid supplied with the index :math:`id=1` (specified in the :ref:`VISUAL_.GRIDS` section of the input) and exporting them to the respective files, use::

  .GRIDFUNCTIONS
  3
  name=ed id_grid=1 purpose=visualization export=cube file_out=ed.cube
  name=ed_laplacian id_grid=1 purpose=visualization export=hdf5 file_out=ed_laplacian.h5
  name=rdg id_grid=1 purpose=visualization export=csv file_out=rdg.csv


Available grid functions
------------------------

The table below summarizes the grid functions currently implemented in `DIRAC` with their respective test directories in `DIRAC` test suite. 


The table below summarizes the grid functions currently implemented in `DIRAC` with their respective test directories. 

name | property | `DIRAC` tests directories
-----|----------|--------------------------
`ed` | electron density, $\rho(\vec{r})$  | `visual_3d_electronden*`
`ed_gradient` | electron density gradient, $\nabla\rho(\vec{r})$ | `visual_3d_electronden*`
`ed_hessian` | electron density Hessian, $H\rho(\vec{r})$   | `visual_3d_electronden*`
`ed_laplacian` | electron density Laplacian, $\nabla^2\rho(\vec{r})$ | `visual_3d_electronden*`
`ed_sign_lambda2` | $sign(\lambda_2)\rho(\vec{r})$ | `visual_3d_electronden*`
`rdg` | reduced density gradient | `visual_3d_electronden*`
`elf` | electron localization function | `visual_3d_elf*`
`kin` | kinetic energy density (1) | `visual_3d_kinden*`
`kinls` | kinetic energy density (2) | `visual_3d_kinden*`
`kinsl` | kinetic energy density (1) | `visual_3d_kinden*`
`kinlap` | the Laplacian of the kinetic energy density | `visual_3d_kinden*`
`kintau` | kinetic energy density (3) | `visual_3d_kinden*`
`kinnr` | the non-relativistic kinetic energy density | `visual_3d_kinden*`
`s` | spin density, $\rho_s(\vec{r})$ | `visual_3d_spinden*`
`divs` | $\Delta \cdot\rho_s(\vec{r})$ | `visual_3d_spinden*`
`rots` | $\Delta \times \rho_s(\vec{r})$ | `visual_3d_spinden*`
`gamma5` | $\gamma^5$ density | `visual_3d_gamma5*`
`edip` | the electric dipole density | `visual_3d_rspE*`
`bdip` | the magnetic dipole density | `visual_3d_rspB*`
`ndip` | the nuclear magnetic dipole density | `visual_3d_rspB*`
`j` | the probability current density, $\vec{j}(\vec{r})$ | `visual_3d_rspB*`
`rotj` | $\nabla \times \vec{j}(\vec{r})$ | `visual_3d_rspB*`
`divj` | $\nabla \cdot \vec{j}(\vec{r})$ | `visual_3d_rspB*`
`gradj` | $\Delta \vec{j}(\vec{r})$ | `visual_3d_rspB*`


Comments:

* the kinetic energy density (1) and its LS- and SL-components (2) are calculated as $c \psi^\dagger(r) (\alpha \cdot p) \psi(r)$
* the kinetic energy density (3) calculated as $\tau = \nabla_i \phi_k \nabla_i \phi_k$



Available `key=value` options
-----------------------------

Notes:

* *current limitations:* 
    * make sure there is one space between each `key=value` entry in the grid function specification line.
    * the maximum line length in the input file is 200 characters

General `key=value` options
~~~~~~~~~~~~~~~~~~~~~~~~~~

* `name=string`

    The name of the grid function to calculate; see the available options in the table above.

* `id_grid=int`

    The ID of a grid on which the grid function should be calculated (see the grid ID specification in :ref:`VISUAL_.GRIDS`.

* `input=option`

    Specify grid function input. Available options:
    
    * `input=calculate` : calculate grid function; the default option
    * `input=import` : import grid function from file

* `purpose=option`

    Specify the type of calculations. Available options:
    
    * `purpose=visualization` : calculate grid function in grid points; possible export as scalar/vector/tensor field for a subsequent visualization; the default option
    * `purpose=integration` : integrate grid function; for the specific examples, see the corresponding test in `DIRAC` test suite: `test/visual_integration`.

* `outpri=option`

    Specify whether the grid function values should also be printed out to `DIRAC` output. Available options: `outpri=no` (default if `purpose=visualization`), `outpri=yes` (default if `purpose=integration`).

* `export=option`

    If the grid function should be exported to a file, this allows to define the file format for export. Supported formats:
    
        * `export=hdf5` : export to hdf5 file
        * `export=txt` : export to text file, no header, space as a separator (for backward compatibility with old VISUAL code)
        * `export=cube` : export to gaussian cube file; available only for 3D rectilinear grids (for backward compatibility with old VISUAL code)
        * `export=csv` : export to text file (csv): header line, coma as a separator

* `file_inp=string`

    The name of a file, from which the grid should be read.

* `file_out=string`

    The name of a file, to which the grid should be exported


`key=value` options allowing for pointwise modification of grid functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* `dscale=f`

    Scale densities *down* by a factor `f`. Default: `dscale=1.0`

* `uscale=f`

    Scale densities *up* by a factor `f`. Default: `uscale=1.0`

* `carpow=[x,y,z]` 

    Scale densities by Cartesian product :math:`x^iy^jz^k`. The `x,y,z` values are three integers specifying the exponents :math:`(i,j,k)`. For example, `carpow=[1,0,0]` is equivalent to calculating the `x`-component of the electric dipole moment density (specification `name=edip`).

* `radpow=f` 

    Scale densities by a radial power :math:`r^{n}`. The keyword is followed by three integers specifying the exponent :math:`n`. Example: `radpow=1` allows to the calculation of radial expectation values :math:`<r>` with respect to the origin.

Notes:

* for the specific examples, see the corresponding test in DIRAC test suite: `test/visual_grid_function_modifications`

Additional keywords affecting grid functions
--------------------------------------------

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

*Warning:* this keyword affects *all* grid functions specified in the input.


keyword(CVALUE)

Set the speed of light to be used in the `VISUAL` module.

*Warning:* this keyword affects *all* grid functions specified in the input.

General setup
-------------

keyword(PRINTL)

Control the print level::

  .PRINTL
  print_level

Predefined `print_level` values:

  * `print_level` = 0: ony the basic information from the `VISUAL` calculations is printed out; the default option
  * `print_level` = `1` : more details on the methods used are printed out
  * `print_level` = `2` : additional print of matrices (warning: large output!)
  * `print_level` > `2` : the development prints (warning: large output!)

All information is printed out to `DIRAC` output.


Additional notes (untested functionalities; notes from the old VISUAL code)
===========================================================================

.. warning::

   Only the calculation of the density is tested for open shell configurations
   (and relies on a correct .OCCUPATION). All other properties are only
   tested for closed shell systems and should not be trusted for open shell systems
   without a thorough testing.



