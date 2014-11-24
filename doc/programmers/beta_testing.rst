

Checklist for beta-testers
==========================


Verify that the tarball installs correctly
------------------------------------------

Check out the release branch::

  $ git checkout release-14

Create the tarball yourself::

  $ ./pamadm --release
  $ cd build
  $ cmake ..
  $ make release

Please report if the tarball generation is broken.

Extract the generated tarball and install it::

  $ ./setup [--fc=ifort --cc=icc --cxx=icpc]

Please report problems. Pay attention whether math libraries are correctly
detected.


Run the test set
----------------

Run the test set using the extracted tarball (see further up)::

  $ ./setup [--flags]
  $ cd build
  $ ctest -jN -D Experimental

The result appears automatically on the
`DIRAC Dashboard <https://repo.theochem.kth.se/cdash/index.php?project=DIRAC>`_.
And we can all see which tests fail and why.


Verify the author list
----------------------

Verify the author list below the logo in the output.  Check versions. Check
typos.


Run tutorials from the website
------------------------------

Verify that they actually run.


Verify the manual and other documentation
-----------------------------------------

Correct typos and mark sections that are wrong
or not clear with a warning::

  .. warning::

     Something wrong here! (replace this text)

This will create a red box with the warning inside.
