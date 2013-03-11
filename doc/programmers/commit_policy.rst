

Commit policy after the release is out
======================================

Since DIRAC12 we have two novelties:
 - We have put a patching/updating mechanism in place for patches and bugfixes.
 - We have a versionned documentation.

Both the patching and the documentation **will** break down if we do not follow
the guidelines outlined here.

1. Commit patches and bugfixes (large or small) that relate to DIRAC12
   to the release-12 branch and not to master.

2. Commit changes to the DIRAC12 documentation to the
   release-12 branch and not to master. In the days
   and weeks following the release we will modify the DIRAC12 documentation a lot.
   For this simply commit to release-12. The documentation
   is generated from that branch. This branch continues to live
   until we release DIRAC13.

To summarize the above two points: all commits that go to master will **never**
make it to this documentation: http://diracprogram.org/doc/release-12/ and will
**never** be part of patches to DIRAC12 if you don't carefully plan your
commits, the patches **will** get lost (from the user's point of view).

**All** commits that go to release-12 can be merged to master at any moment.
This means that all changes to release-12 will get integrated into the main line
development.
**Never** merge from master to release-12.

If, by accident, you committed to the wrong branch, you can ``git cherry-pick`` the
commit to the release-12 branch.

3. Do not change the file VERSION. This file is changed when
   we release a new patch. But we will not release a patch from each
   commit to release-12

4. If you commit a bugfix or patch that is important please send Radovan an email
   so that he can bump the VERSION and publish a new patch.

To summarize, ask yourself the following questions:
 - Is this code or documentation change useful for DIRAC12 users? If yes, commit to release-12.
 - Is this code or documentation change only useful for future versions of DIRAC (13 or later)? If yes, commit to some other branch.
