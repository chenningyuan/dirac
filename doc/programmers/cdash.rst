

Nightly testing dashboard
-------------------------

Before you start, set::

  --scratch=/your/scratch/path

in ``~/.diracrc``.
In the case of a parallel DIRAC nightly run, you may have to set up more
parameters, like::

  --machfile=$PBS_NODEFILE
  --global-scratch-disk

It is very easy to contribute to the
`DIRAC Dashboard <http://repo.ctcc.no/CDash/index.php?project=DIRAC>`_.
You don't need to install anything extra except CMake.
It is as easy as typing::

  $ ./setup [--flags] -D BUILDNAME="my-build-name"
  $ cd build
  $ make Experimental

The result of ``make Experimental`` appears automatically on the
`DIRAC Dashboard <http://repo.ctcc.no/CDash/index.php?project=DIRAC>`_.

Once ``make Experimental`` works you can create a script to run ``make Nightly``.
To get inspired have a look at the cdash example recipes in the repository
under ``maintenance/cdash``. Then put it under crontab. Here is an example::

  # m  h  dom mon dow command
    05 00 *   *   *   /home/bast/nightly/smokey_debian-6.0_gfortran-4.4.5
    00 02 *   *   *   killall -9 dirac.x
    05 02 *   *   *   /home/bast/nightly/smokey_debian-6.0_gfortran-4.4.5_i8
    00 04 *   *   *   killall -9 dirac.x
    05 04 *   *   *   /home/bast/nightly/smokey_debian-6.0_gfortran-4.4.5_i8_mpi
    00 06 *   *   *   killall -9 dirac.x
