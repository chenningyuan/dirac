
Contributing to the DIRAC Dashboard
-----------------------------------

Before you start, set --scratch=/your/scratch/path in ~/.diracrc.

In the case of parallel DIRAC installation, you have to set up more
parameters in the diracrc (=~/.diracrc, or .diracrc or diracrc) file, like
--machfile=$PBS_NODEFILE
--scratch=/your/scratch/path
--global-scratch-disk

It is very easy to contribute to the DIRAC Dashboard:
http://repo.ctcc.no/CDash/index.php?project=DIRAC

You don't need to install anything extra except CMake.
It is as easy as typing::

  $ ./setup [--flags] -D BUILDNAME="my-build-name"
  $ cd build
  $ make Experimental

The result of "make Experimental" appears automatically on the Dashboard.

Make Experimental takes a long time (2 hours or more). This is because we have
many very slow tests. For experimenting with make experimental you can cheat
and remove many lines in cmake/Tests.cmake (and do not rerun pamadm).  Later
you can regenerate the full cmake/Tests.cmake.

Once "make Experimental" works you can create a script to run "make Nightly".
To get inspired have a look at the cdash example recipes in the repository
under maintenance/cdash. Then put it under crontab. Here is an example::

  # m  h  dom mon dow command
    05 02 *   *   *   /home/bast/nightly/smokey_debian-6.0_gfortran-4.4.5
    00 04 *   *   *   killall -9 dirac.x
    05 04 *   *   *   /home/bast/nightly/smokey_debian-6.0_gfortran-4.4.5_i8
    00 06 *   *   *   killall -9 dirac.x
    05 06 *   *   *   /home/bast/nightly/smokey_debian-6.0_gfortran-4.4.5_i8_mpi
    00 08 *   *   *   killall -9 dirac.x
