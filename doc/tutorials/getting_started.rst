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

Dirac-Coulomb SCF
-----------------

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

Closed shell CCSD
-----------------

 Running CCSD(T) calculations is straightforward if the ground state of
the molecule can be well-described by a single closed shell determinant.
This is the case for many molecules for which the input requires
specification of the basis set (TZ or better is recommended, but for
this example we will take DZ to reduce the run time). We take the inter
halogen molecule ClF as an example and use the default cut-offs for
correlating electrons (include all valence electrons with energy above
-10 Hartree) and truncation of virtual space (delete virtuals above 20
Hartree). The geometry (**clf.xyz**) file is very simple, requiring only
specification of the coordinates. The program will then identity the
symmetry as |C\_{\\inf v}|:

::

    2
    ClF molecule at equilibrium distance taken from NIST.
    Cl   0.0                 0.0                 0.0
    F    0.0                 0.0                 1.628

We specify the wave function type (SCF, followed by RELCCSD) and basis
set in the input file (**cc.inp**) to calculate the CCSD(T) energy with
the Dirac-Coulomb Hamiltonian (with the contribution from SS\|SS
integrals approximated by a simple Coulombic correction).

::

    **DIRAC
    .WAVE FUNCTION
    **WAVE FUNCTION
    .SCF
    .RELCCSD
    **MOLECULE
    *BASIS
    .DEFAULT
    cc-pVDZ
    *END OF INPUT

 Now start the calculation by executing:

::

    pam  --mol=clf.xyz --inp=cc.inp

If everything works fine the results of the calculation will be written
to **cc\_clf.out**. In addition an archive file *cc\_clf.tgz* will be
created, as mentioned above. 

You can suppress creation of the archive file by pam flag  *--noarch*.
