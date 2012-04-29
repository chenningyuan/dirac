Closed shell CCSD
=================

From Diracwiki
~~~~~~~~~~~~~~

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

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
created. Depending on the type of calculation the archive file may
include:

-  ``DFCOEF`` containing molecular orbital coefficients and possibly
   energies.
-  ``*.cube`` containing plots of molecular properties, typically
   electron densities, in the gridded Gaussian Cube format.
-  Additional files for restarting particular types of calculations.

