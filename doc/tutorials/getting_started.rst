Getting started with DIRAC
==========================

DIRAC uses two files as inputs: the input file **\*.inp**, which
determines what should be calculated, and either a geometry file
**\*.xyz**, that provides the nuclear coordinates, or a molecule file
**\*.mol**, which defines the molecular geometry and also the basis set.
These files are described in detail in the
`Manual </dirac/index.php/Manual>`_. 

When DIRAC runs it writes its
output to a text file, named after the molecule and input file names. If
this file already exists, for example from a previous calculation, the
old file is renamed to make place for the new output. The **pam** python
script that launches the main DIRAC executable also writes some general
information to standard output.

Let's try it out and run our first DIRAC calculation with the following
minimal input (**hf.inp**) to calculate the Dirac-Coulomb Hartree-Fock
ground state (with the default simple Coulombic correction
(`LVCORR </dirac/index.php/Manual:HAMILTONIAN#.LVCORR>`_) to approximate
the contribution from SS\|SS integrals, see [Visscher1997a]_ ::

  **DIRAC
  .WAVE FUNCTION
  **WAVE FUNCTION
  .SCF
  **MOLECULE
  *BASIS
  .DEFAULT
  cc-pVDZ
  *END OF INPUT

together with the following example geometry file (**methanol.xyz**)::

  6
  Methanol. This is a comment line.
  C       0.000000000000       0.138569980000       0.355570700000   
  O       0.000000000000       0.187935770000      -1.074466460000  
  H       0.882876920000      -0.383123830000       0.697839450000  
  H      -0.882876940000      -0.383123830000       0.697839450000  
  H       0.000000000000       1.145042790000       0.750208830000  
  H       0.000000000000      -0.705300580000      -1.426986340000

Now start the calculation by executing::

  pam  --mol=methanol.xyz --inp=hf.inp

If everything works fine the results of the calculation will be written
to **hf\_methanol.out**. In addition an archive file *hf\_methanol.tgz*
will be created. Depending on the type of calculation the archive file
may include:

-  ``DFCOEF`` containing molecular orbital coefficients and possibly
   energies.
-  ``*.cube`` containing plots of molecular properties, typically
   electron densities, in the gridded Gaussian Cube format.
-  Additional files for restarting particular types of calculations.
