Getting started with DIRAC
==========================

You have installed and tested DIRAC and now it's time to run our first DIRAC
calculation!

You have two possibilities to run DIRAC. The traditional uses two input files:
the inp file, which determines what should be calculated, together with the mol
file, which defines the molecular geometry and the basis set.  Alternatively
you can use the inp file together with a geometry file (xyz file). In this case
the xyz provides the nuclear coordinates and the inp file contains in addition
to the job description also the basis set information (see :ref:`molecule`).

DIRAC uses two files as inputs: the input file **\*.inp**, which
determines what should be calculated, and the molecule file **\*.mol**
or **\*.xyz** , which defines the molecular geometry and possibly the
basis set. These files are described in detail in the
`Manual </dirac/index.php/Manual>`_. When DIRAC runs it writes its
output to a text file, named after the molecule and input file names. If
this file already exists, for example from a previous calculation, the
old file is renamed to make place for the new output. The **pam** script
that launches the main DIRAC executable also writes some general
information to standard output.

DIRAC writes its output to a text file, named after the molecule and input file
names. If this file already exists, for example from a previous calculation,
the old file is renamed to make place for the new output. The python script
pam that launches the main DIRAC executable also writes some general
information to standard output.

Let's try it out and run our first DIRAC calculation with the following
minimal input (**dc\_hf.inp**) to calculate the Dirac-Coulomb
Hartree-Fock ground state

::

    **DIRAC
    .WAVE FUNCTION
    **WAVE FUNCTION
    .SCF
    *END OF INPUT

together with the following example molecule file (**ch3f.mol**):

::

    DIRAC
    CH3F - this line is a comment
    this line is also a comment 
    C   3              
            6.    1
    C    0.0                 0.0                 0.426596314746
    LARGE BASIS cc-pVDZ
            9.    1
    F    0.0                 0.0                -0.98392891998
    LARGE BASIS cc-pVDZ
            1.    3
    H1   0.509715617414      0.882853346772      0.786416057331
    H2   0.509715617414     -0.882853346772      0.786416057331
    H3  -1.01943123483       0.0                 0.786416057331
    LARGE BASIS cc-pVDZ
    FINISH

Now start the calculation by executing:

::

    pam --mb=200 --mol=ch3f.mol --inp=dc_hf.inp

The **--mb=200** flag specifies that DIRAC should allocate 200 MB for
the work array. If everything works fine the results of the calculation
will be written to **dc\_hf\_ch3f.out**. In addition an archive file
**dc\_hf\_ch3f.tgz** will be created. Depending on the type of
calculation the archive file may include:

-  ``DFCOEF`` containing molecular orbital coefficients and possibly
   energies.
-  ``*.cube`` containing plots of molecular properties, typically
   electron densities, in the gridded Gaussian Cube format.
-  Additional files for restarting particular types of calculations.


Dirac-Coulomb SCF
-----------------

Let's try it out and run our first DIRAC calculation with the following minimal
input (hf.inp) to calculate the Dirac-Coulomb Hartree-Fock ground state
(with the default `simple Coulombic correction
<http://dx.doi.org/10.1007/s002140050280>`_ to approximate the contribution
from (SS|SS) integrals::

  **DIRAC
  .WAVE FUNCTION
  **WAVE FUNCTION
  .SCF
  **MOLECULE
  *BASIS
  .DEFAULT
   cc-pVDZ
  *END OF INPUT

together with the following example geometry file (methanol.xyz)::

  6
  my first DIRAC calculation # anything can be in this line
  C       0.000000000000       0.138569980000       0.355570700000   
  O       0.000000000000       0.187935770000      -1.074466460000  
  H       0.882876920000      -0.383123830000       0.697839450000  
  H      -0.882876940000      -0.383123830000       0.697839450000  
  H       0.000000000000       1.145042790000       0.750208830000  
  H       0.000000000000      -0.705300580000      -1.426986340000

Now start the calculation by executing::

  pam --mol=methanol.xyz --inp=hf.inp

If everything works fine the results of the calculation will be written
to hf\_methanol.out. In addition an archive file hf\_methanol.tgz
will be created. Depending on the type of calculation the archive file
may include:

-  DFCOEF containing molecular orbital coefficients and possibly
   energies.
-  cube files containing plots of molecular properties, typically
   electron densities, in the gridded Gaussian Cube format.
-  Additional files for restarting particular types of calculations.


Closed shell CCSD
-----------------

Running CCSD(T) calculations is straightforward if the ground state of
the molecule can be well-described by a single closed shell determinant.
This is the case for many molecules.

The input requires the
specification of a basis set (TZ or better is recommended, but for
this example we will take DZ to reduce the run time). We take the inter
halogen molecule ClF as an example and use the default cut-offs for
correlating electrons (include all valence electrons with energy above
-10 hartree) and truncation of virtual space (delete virtuals above 20
hartree). The geometry file (clf.xyz) is very simple and only requires
to specify the coordinates. The program will then identify the
symmetry as :math:`C_{\infty v}`::

  2
  ClF molecule at equilibrium distance taken from NIST.
  Cl   0.0  0.0  0.0
  F    0.0  0.0  1.628

We specify the wave function type (SCF, followed by RELCCSD) and basis
set in the input file (cc.inp) to calculate the CCSD(T) energy with
the Dirac-Coulomb Hamiltonian again with the contribution from (SS|SS)
integrals approximated by a simple Coulombic correction::

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

Now start the calculation by executing::

  pam --mol=clf.xyz --inp=cc.inp

If everything works fine the results of the calculation will be written
to cc_clf.out. In addition an archive file cc_clf.tgz will be
created, as mentioned above. 

You can suppress creation of the archive file by::

  pam --noarch
