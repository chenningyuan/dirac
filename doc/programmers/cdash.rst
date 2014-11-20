

Nightly testing dashboard
-------------------------

It is very easy to contribute to the
`DIRAC Dashboard <https://repo.theochem.kth.se/cdash/index.php?project=DIRAC>`_.
You don't need to install anything extra except CMake.
It is as easy as typing::

  $ ./setup [--flags]
  $ cd build
  $ make -jN -D Nightly [--track master|--track release-14]

The result of appears automatically on the
`DIRAC Dashboard <https://repo.theochem.kth.se/cdash/index.php?project=DIRAC>`_.

To get inspired have a look at the script cdash_crontab.sh recipes in the repository
under ``maintenance/cdash``. Then put it under crontab. Here is an example::

  # m  h  dom mon dow command
    01 03 *   *   *   /home/user/nightly/dirac/maintenance/cdash/cdash_crontab.sh -b 'master'     -p '--adm'     -m 12      -c 'Intel' -n 'Intel-15.0.0'         -t 'master'
    20 03 *   *   *   killall -9 dirac.x
    21 03 *   *   *   /home/user/nightly/dirac/maintenance/cdash/cdash_crontab.sh -b 'release-14' -p '--release' -m 12 -r 1 -c 'Intel' -n 'Intel-15.0.0'         -t 'release-14'
    40 03 *   *   *   killall -9 dirac.x
    41 03 *   *   *   /home/user/nightly/dirac/maintenance/cdash/cdash_crontab.sh -b 'release-14' -p '--release' -m 12 -r 1 -c 'Intel' -n 'Intel-15.0.0-i8' -i 1 -t 'release-14'
    59 03 *   *   *   killall -9 dirac.x
    01 04 *   *   *   /home/user/nightly/dirac/maintenance/cdash/cdash_crontab.sh -b 'master'     -p '--adm'     -m 12      -c 'GNU'   -n 'GNU-4.8.2'         -t 'master'
    20 04 *   *   *   killall -9 dirac.x
    21 04 *   *   *   /home/user/nightly/dirac/maintenance/cdash/cdash_crontab.sh -b 'release-14' -p '--release' -m 12 -r 1 -c 'GNU'   -n 'GNU-4.8.2'         -t 'release-14'
    40 04 *   *   *   killall -9 dirac.x
    41 04 *   *   *   /home/user/nightly/dirac/maintenance/cdash/cdash_crontab.sh -b 'release-14' -p '--release' -m 12 -r 1 -c 'GNU'   -n 'GNU-4.8.2-i8' -i 1 -t 'release-14'
    59 04 *   *   *   killall -9 dirac.x
