

The scratch directory is the place where DIRAC will write temporary files. At
the beginning of the run, DIRAC copies input files and the binary into the
scratch space and at the end of the calculation, the output is copied back,
possibly together with other useful (restart) files.


How you can set the scratch directory
=====================================

You can set the scratch space like this::

  $ pam --scratch=/somepath

If ``--scratch`` is not set, DIRAC will look for the environment variabls
``DIRAC_TMPDIR``.

If these variables are not set, DIRAC defaults to
``$HOME/DIRAC_scratch_directory``.

On your local computer you typically want to use the same scratch space every
time and you may be tired of typing ``pam --scratch`` every time.  In this
case, put ``--scratch=/somepath`` in your ``~/.diracrc``.

DIRAC will always append ``$USER/DIRAC_inp_mol_pid`` to the scratch path.  This
is to prevent users from accidentally wiping out their home directory.


What if you don't want DIRAC to append anything to your scratch path
====================================================================

For this use::

  $ pam --scratchfull=/somepath

In this case DIRAC will use the path as it is and not append anything to it.
Be careful with this since DIRAC may remove the path after the calculation
has finished.


How to set the scratch path on a cluster
========================================

On a cluster you typically want to be able to find the scratch directory based
on some job id. If you are on a cluster with a Torque/PBS scheduler then it can
be useful to set::

  $ pam --scratch=/somepath/$USER/$PBS_JOBID


Cluster with a global scratch disk
----------------------------------

This means that all nodes write temporary files to the same disk.

No additional flags are needed for this run.


Cluster with local scratch disks
--------------------------------

This means that each node writes temporary files to its own disk.

For this situation the input files have to be distributed over
the nodes based on ``--machfile``::

  $ pam --scratch=/somepath/$USER/$PBS_JOBID --machfile=/path/to/machfile

If you don't give ``--machfile``, DIRAC assumes that you run using a global
scratch disk.

The list of machines is typically generated when the scheduler launches
the job but you can also set the list of machines manually.

The input files are by default distributed using rsh/rcp. If your cluster
does not support these you can change to ssh/scp protocol::

  $ pam --scratch=/somepath/$USER/$PBS_JOBID --machfile=/path/to/machfile --rsh=ssh --rcp=scp

Replace strings in input files
==============================

The " pam --replace" 
is a fantastic option (i can say this because
it's Ulf's idea) and I encourage you to try it
(before the wrapper.py starts to sue pam for patent infringement).

Try this:

::
 **DIRAC
 .WAVE FUNCTION
 **HAMILTONIAN
 .hamiltonian
 **WAVE FUNCTION
 .SCF
 *END OF

and this: 

::
 DIRAC
 .
 .
 C 1
 10. 1
 Ne 0.0 0.0 0.0
 LARGE BASIS basis
 FINISH

now run it with:

::
 pam --replace hamiltonian=LEVY-LEBLOND --replace basis=cc-pVDZ

it will write to the output file:

::
 inp_mol_hamiltonian=LEVY-LEBLOND_basis=cc-pVDZ.out

you get the point: you can now run whole tables with just two input files!

In another example, to make the mol- input file completely modifiable (generic), 
one replaces the atom name, the proton number and the basis set name at once:

::
 DIRAC 
 . 
 . 
 C 1 
 proton_nu 1 
 atom_name 0.0 0.0 0.0 
 LARGE BASIS basis 
 FINISH 

This launching script command shall be:

::
 pam --replace proton_nu=10.0  --replace atom_name=Ne --replace basis=cc-pVDZ 

A nice application of this is to make a PES scan, written as a script:

::
 for r in `seq 1.00 0.05 2.00`; do
    pam --replace distance=$r ..
 done
