

Testing the installation
========================

It is very important that you verify that your DIRAC installation correctly
reproduces the reference test set before running any production calculations.

The test set driver is CTest which can be invoked with "make test" after building
the code.


Environment variables for testing
---------------------------------

Before testing with "make test" you should export the
following environment variables::

  $ export DIRAC_TMPDIR=/scratch        # scratch space (adapt the path of course)
  $ export DIRAC_NUM_MPI_PROCS=4        # in this case 4 processes, only relevant if you compile with MPI

Note that if you set the DIRAC_NUM_MPI_PROCS to something different from 1,
the pam script will assume you have compiled using MPI and run the mpirun
command!


Running the test set
--------------------

You can run the whole test set either using::

  $ make test

or directly through CTest::

  $ ctest

Both are equivalent ("make test" runs CTest) but running
CTest directly makes it easier to run sequential tests on several
cores::

  $ ctest -j4        # only for sequential tests

You can select the subset of tests by matching test names to a regular expression::

  $ ctest -R dft

or matching a label::

  $ ctest -L quick

To print all labels::

  $ ctest --print-labels
