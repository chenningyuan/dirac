

How to define an alternative MPI launcher
-----------------------------------------

The pam script by default launches MPI runs with mpirun. If you cannot use mpirun
because you are on a system which does not support it or which offers alternative MPI launchers
(SGI, Cray), you can define your own DIRAC_MPI_COMMAND like this::

  $ export DIRAC_MPI_COMMAND="mpirun -np 8"
  $ export DIRAC_MPI_COMMAND="mpprun"
  $ export DIRAC_MPI_COMMAND="aprun -n 24"

This will then launch::

  $DIRAC_MPI_COMMAND dirac.x

instead of::

  mpirun -np N dirac.x
