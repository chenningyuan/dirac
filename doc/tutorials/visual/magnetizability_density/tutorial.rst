

Magnetizability density
=======================

In this tutorial we will plot the CGO magnetizability density of the water
molecule.  We will do this in two steps: first we will run the SCF and response
calculation and save the DFCOEF and PAMXVC files, then we will restart from
these and calculate the magnetizability density on a grid of points and write
this into a cube file which can be opened by your favorite molecular
visualization program (the author of this tutorial likes Avogadro).


Getting the wave function and the response vectors
--------------------------------------------------

Let us use the following molecule input (h2o.mol):

.. literalinclude:: h2o.mol

Together with a simple job input (magnetizability.inp):

.. literalinclude:: magetizability.inp

We can get the wave function with the following run script:

.. literalinclude:: magetizability.run

After running the script you should see the files DFCOEF and PAMXVC in your
submit directory.


Plotting the magnetizability density
------------------------------------

We will use the following run script to generate the cube file:

.. literalinclude:: density.run

Together with the following job input (density.inp):

.. literalinclude:: density.inp

In the output we should verify that we integrate to the correct
isotropic magnetizability::

     scalar              x-component         y-component         z-component
    
  -0.1710573104E+01    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00

Finally you can open plot.3d.cube with your favorite program and visualize the
density. By selectively deactivating contributions in density.inp you can
single out specific components or paramagnetic and diamagnetic contributions
separately.
