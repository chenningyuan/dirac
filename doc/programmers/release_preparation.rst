




Release preparation
===================


Feature freeze
--------------

Two months before the release (example: DIRAC12), a release branch (example:
release-12) is created from master.  This is a feature freeze and from this
moment on the release branch only receives cosmetics and bugfixes, no new
features, no merges from master or other branches (except cherry-picks; more
about it later).  You as developer are responsible to commit or merge the
to-be-released code to master, before the release branch is created.  If your
code is part of a pamadm module and you want to release it, make sure that it
is part of pamadm --release. If you do not want to release your pamadm module,
make sure it is not part of pamadm --release.


Can I commit bleeding edge code to master after the feature freeze?
-------------------------------------------------------------------

Yes you can! Exactly for this reason we have the release branch and exactly for
this reason we will never merge from master to the release branch, only from
the release branch to master.


How to check out the release branch (example release-12 for DIRAC12)
--------------------------------------------------------------------

Check out the release branch::

  $ git checkout release-12

Now you have it next to master, try it out::

  $ git branch

You can switch back to master::

  $ git checkout master

And then switch back to the release branch::

  $ git checkout release-12


How to commit changes and bugfixes that affect the soon-to-be-released code
---------------------------------------------------------------------------

Commit all such changes to the release branch, not to master.  This way also no
commit or bugfix will get lost.  Please do *not* transfer commits from the
release branch to master manually.  This is not only unnecessary work but it is
harmful (conflicts)!
Again, do not commit to master, commit to the release branch::

  $ git checkout release-12
  $ git pull                      # update the local release-12 branch
  $ git commit                    # commit your modifications
  $ git push                      # push changes to release-12
  $ git checkout master
  $ git pull                      # update your local master branch
  $ git merge --no-ff release-12  # merge changes to master
  $ git push                      # push changes to master

If you have no time, you don't have to merge to master yourself, you can wait
until somebody else does it for his/her commits. No commit can get lost.

O.K. what if you accidentally committed something to master which belongs to
the release? In this case do *not* merge master to the release branch, but
rather cherry-pick the individual commit(s) to the release branch::

  $ git checkout release-12
  $ git cherry-pick [hash]        # with [hash] that corresponds
                                  # to the commit


How to exclude code from being released
---------------------------------------

You as the developer of the are responsible for removing your code from the
release branch that you do not wish to be released.  This is *not* done by
actually removing the code but by using pamadm in combination with CPP
statements. So again, *do not actually delete* code lines, use CPP statements
as demonstrated below:

Modules can be excluded from the release version by surrounding the appropriate
parts of code by preprocessor statements, e.g.::

  #ifdef MOD_KRMC
           CALL TKRMCORB
           CALL SETDC2(0)
           CALL PSIOPT('MCSCF',WFCONV,WORK,LWORK)
           IF(DOANA) CALL PAMANA(WORK,LWORK)
  #else
              CALL QUIT(
       &  'Second order KMCSCF optimization not included in this version')
  #endif

This means that the calls are included only if MOD_KRMC is defined.  Look in
pamadm for existing modules and how they define their CPP flags and
directories.  If your module is not part of pamadm --release, then the code
will be automatically removed by the release script. In the above example, the
script would remove the actual MOD_KRMC code and what would be left is::

              CALL QUIT(
       &  'Second order KMCSCF optimization not included in this version')

If you want to protect code from being released which is not part of any pamadm
module, you can use MOD_UNRELEASED::

  #ifdef MOD_UNRELEASED
  !      secret code
  #endif

To summarize, you can remove code using either existing CPP filters or using
MOD_UNRELEASED. Do this on the release branch and merge to master::

  $ git checkout release-12
  $ git pull                      # update the local release-12 branch
  $ git commit                    # commit your modifications
  $ git push                      # push changes to release-12
  $ git checkout master
  $ git pull                      # update your local master branch
  $ git merge --no-ff release-12  # merge changes to master
  $ git push                      # push changes to master


How to release or export the code
---------------------------------

Most important rule: treat unpublished code of others with respect.  Do not
release code of others without asking.

1. Verify copyright headers.
2. Verify logo (list of authors).
3. Bump the version.

As a final step create the tarball::

  $ ./pamadm --release
  $ mkdir build
  $ cd build
  $ cmake ..
  $ make release
