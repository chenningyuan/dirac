

Nightly testing dashboard
-------------------------

Please contribute to the
`DIRAC Dashboard <https://testboard.org/cdash/index.php?project=DIRAC>`__.
You don't need to install anything extra except CMake.
It is as easy as typing::

  $ ./setup [--flags]
  $ cd build
  $ make -jN -D Nightly [--track master|--track release-14]

The result of appears automatically on the
`DIRAC Dashboard <https://testboard.org/cdash/index.php?project=DIRAC>`__.

To get inspired have a look at the script ``cdash_crontab.sh`` in the repository
under ``maintenance/cdash`` - this is the script
that Radovan uses. Then put it under crontab. Here is an example::

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

We currently have the following tracks installed: master, release-14, Nightly,
Experimental, BoundsCheck, Grid, Miro, and Coverage.  It is no problem to
install additional tracks - just ask Radovan.  The tracks are there to organize
the various builds. So instead of submitting a build to the track "Nightly"
with the build name "master-foo" it is clearer and more concise to submit a
build "foo" to the track "master".  The tracks "BoundsCheck", "Grid", and
"Miro" are used by Miro.  I (Radovan) recommend to use "master" or "release-14"
depending on the context and not to use "Nightly" or "Experimental". But it is
OK to use "Nightly" or "Experimental" to calibrate your scripts.

The semantic difference between "Nightly" and "Experimental" is that "Nightly"
corresponds to a Git hash closest to a certain configurable time whereas
"Experimental" corresponds to now. This means that for extremely actively
developed projects it may be important to test all nightly tests using a well
defined version and then one should prefer "Nightly" over "Experimental".
