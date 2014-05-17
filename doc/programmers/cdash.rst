

Nightly testing dashboard
-------------------------

It is very easy to contribute to the
`DIRAC Dashboard <https://repo.theochem.kth.se/cdash/index.php?project=DIRAC>`_.
You don't need to install anything extra except CMake.
It is as easy as typing::

  $ ./setup [--flags]
  $ cd build
  $ make -jN -D Nightly

The result of appears automatically on the
`DIRAC Dashboard <https://repo.theochem.kth.se/cdash/index.php?project=DIRAC>`_.

To get inspired have a look at the cdash example recipes in the repository
under ``maintenance/cdash``. Then put it under crontab. Here is an example::

  # m  h  dom mon dow   command
    05 01 *   *   *     /home/bast/nightly/dirac/maintenance/cdash/Linux-x86_64-GNU-SYSTEM_NATIVE         -b release-13 -p "--release" -n 12
    00 02 *   *   *     killall -9 dirac.x
    05 02 *   *   *     /home/bast/nightly/dirac/maintenance/cdash/Linux-x86_64-Intel-MKL                 -b release-13 -p "--release" -n 12
    00 03 *   *   *     killall -9 dirac.x
    05 03 *   *   *     /home/bast/nightly/dirac/maintenance/cdash/Linux-x86_64-OpenMPI-GNU-SYSTEM_NATIVE -b release-13 -p "--release" -n 6
    00 05 *   *   *     killall -9 dirac.x
    05 05 *   *   *     /home/bast/nightly/dirac/maintenance/cdash/Linux-x86_64-OpenMPI-Intel-MKL         -b release-13 -p "--release" -n 6
    00 07 *   *   *     killall -9 dirac.x
