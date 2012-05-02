Tutorial:Getting started
========================

From Diracwiki
~~~~~~~~~~~~~~

(Redirected from `Manual:Running
Dirac </dirac/index.php?title=Manual:Running_Dirac&redirect=no>`_)

Jump to: `navigation <#column-one>`_, `search <#searchInput>`_

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

