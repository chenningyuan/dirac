Restarting SCF calculations
===========================

A DIRAC calculation can be restarted at various stages depending on the type of
calculation. In general it should be enough to give the archive file of a
previous calculation as an argument to the ``pam`` script.  DIRAC will then
find and use the available data for restart.  Alternatively you can copy
restart files in using ``pam --put`` and if DIRAC finds them in the scratch
directory, it will automatically restart from them.  In general the program can
restart from

-  ``DFCOEF`` Contains orbital coefficients in machine-readable format, and will be used for
   restarting the SCF procedure.
-  ``DFPCMO`` Contains orbital coefficients in human-readable format, and will be used for
   restarting the SCF procedure.
-  ``PAMFCK`` The Fock matrix. Will be used for
   restarting the SCF procedure.
-  ``AOPROPER`` Contains one-electron integrals.
-  ``X2CMAT`` Contains four- to two-component transformation matrices,
   and allows you to skip the transformation step of a two component
   calculation.


Restarting X2C calculations
===========================

Here we will show how one can run a Hartree-Fock (HF) calculation in
2-component mode with the exact two component (X2C) Hamiltonian and
restart a property calculation at the HF level in a second run without
the need to re-do the X2C transformation.

For the initial SCF step we use the following input file (copy the content below into the ``scf_x2c.inp`` file)::

  **DIRAC
  .WAVE FUNCTION
  **MOLECULE
  *BASIS
  ! Use the same basis for all atoms
  .DEFAULT
  dyall.av3z
  **HAMILTONIAN
  ! Use the exact two-component Hamiltonian
  ! with atomic mean field spin-orbit
  .X2C
  **INTEGRALS
  *READIN
  ! A limitation in the AMFI code requires an uncontracted basis
  .UNCONTRACT
  **WAVE FUNCTION
  .SCF
  *END OF

together with the following molecule file (copy it into the ``br2.xyz`` file)::

  2
  Br 0.0 0.0  1.244
  Br 0.0 0.0 -1.244

This molecule has linear symmetry (:math:`D_{\infty h}`),
and the program will detect this automatically. In most parts of the
program the highest abelian subgroup (:math:`D_{2 h}`) will be used instead of the
full linear symmetry. In this example we use Dyall's augmented TZ basis
dyall.av3z in uncontracted form.

Prior to the HF calculation, the Hamiltonian transformation from 4-
to 2-components will be performed and for restart purposes one should
therefore save the files ``X2CMAT`` and ``AOMOMAT`` containing the
Hamiltonian information and the file ``DFCOEF`` containing the
(converged) MO-coefficients.

We run the initial SCF step using::

  pam --inp=scf_x2c.inp --mol=br2.xyz --get="X2CMAT AOMOMAT DFCOEF" 

The output will be written to the file ``scf_x2c_br2.out``, which can be
monitored during the calculation.

In order to proceed in a second step with a property calculation of,
e.g. the contact density at the Br nucleus at the HF level, we use the following
input file (copy and paste it into the ``property_x2c.inp`` input file)::

  **DIRAC
  .PROPERTIES
  **MOLECULE
  *BASIS
  ! Use the same basis for all atoms
  .DEFAULT
  dyall.av3z
  **HAMILTONIAN
  ! Use the exact two-component Hamiltonian
  ! with atomic mean field spin-orbit
  .X2C
  **INTEGRALS
  *READIN
  ! A limitation in the AMFI code requires an uncontracted basis
  .UNCONTRACT
  **PROPERTIES
  .RHONUC
  *END OF

Restarting from the
(converged) MO-coefficients and reading the transformed X2C Hamiltonian
from file is done automatically once the files are present. We therefore
restart the calculation using::

  pam --inp=property_x2c.inp --mol=br2.xyz --put="X2CMAT AOMOMAT DFCOEF" 

How to run a Coupled-Cluster calculations in separate steps
===========================================================

Most of the time there are no problems in executing the three steps in a
coupled-cluster calculation with RELCCSD: 1) the Hartree-Fock calculation, 2)
the MOLTRA step, the integral transformation from atomic to molecular orbitals, and 3) the
coupled-cluster calculations itself.

However, in some cases it is not the best approach to perform all these in one
go. This is especially true for calculations that use a lot of memory in the
coupled-cluster part but are more economical in the previous steps (say, if the
system has little or no symmetry or, in Fock-Space calculations, if the active
space is large).

In such cases, it is possible to do one step, save the results necessary
for the next calculation, and then run the next step. This can
considerably shorten the total wall time spent on a parallel
calculation, since a larger number of MPI processes can be started for
the SCF (and/or MOLTRA) step(s) with the same amount of memory per
compute node.


Step 1: SCF
-----------

From the SCF part, one has to save the ``DFCOEF`` file after a successful
calculation. If using the 2-component formalism, one can also
save the files ``X2CMAT`` and ``AOMOMAT`` (and possibly ``AOPROPER``), in order to
skip the generation of the 4-component to 2-component transformation
matrices. If your system supports it, you could keep the scratch
directory and then avoid the copying back and forth of these files.

Currently DIRAC's execution script automatically retrieves these
files and places them in a gzip-compressed tar archive, so it is probably not
necessary to retrieve the files individually if this archive is kept.


Step 2: MOLTRA
--------------

MOLTRA will need the ``DFCOEF`` file, and also the one-electron Fock matrix (if
the 2-component methods are used and also in the case of an frozen density
embedding calculation, this will have to be re-generated, so that's why it is a
good idea to keep the ``X2CMAT`` and ``AOMOMAT`` files around.

The result of the MOLTRA step will be a set of files that start with ``MDC...``.
For **SCHEME 4** these files can be considered as 'independent' of
the number of processes used in the parallel run (because only one,
created by the master, will have information used by other modules like
RELCCSD - such as the number of spinors transformed etc.), so it does not
matter how the calculation from which they were obtained was performed.

This is not true for **SCHEME 6**, because there the information
is needed by other modules, and then one has to be careful to use the same
number of MPI processes in the MOLTRA step and beyond.


Step 3: RELCCSD
---------------

RELCCSD will need the ``MDC*`` files generated by MOLTRA. One must be
careful in the input, however; the keyword ``.NO4INDEX`` has to be
specified, and the ``**MOLTRA`` section must not be present in the input.
Again, for **SCHEME 6** one must be careful to have the same number of MPI processes as
in the MOLTRA step.


A note on Fock-Space Calculations
---------------------------------

In the case of Fock-Space calculations, an additional trick can be used: if
previous/lower sectors are converged, RELCCSD can be restarted and be made skip
them. In order to do that one has to do the following apart from the actions
outlined above: First, instead of ``MAXIT=n``, use ``MAXITk=n`` (where k
denotes the sector: 00, 10, 01, 11, 02, 20) and set to zero all the ``MAXITab``
for converged sectors to zero. Then restart exactly the same calculation as
before (in terms of number of nodes etc.). It is important to have the
intermediate files for RELCCSD present in the scratch directory in this case.
The easiest way to do so is to keep the scratch directory from the previous
run.


Troubleshooting
===============

Can I move the DFCOEF file from one machine to another and restart from it?
---------------------------------------------------------------------------
                                                                           
In general no, unless it is the same architecture.

What you can do instead is to export the DFCOEF in a machine independent
format using the .PCMOUT keyword. This creates a file called DFPCMO.
Copy DFPCMO to the other machine, copy it to the scratch
directory using the pam script, and DIRAC will correctly restart from it.


I tried to restart but get the error "Incompatible number of basis functions" - what does it mean?
--------------------------------------------------------------------------------------------------

This probably means that you have changed the basis set. Presently DIRAC cannot
restart from DFCOEF or DFPCMO if the basis set changes (more correctly if
matrix dimensions change).
