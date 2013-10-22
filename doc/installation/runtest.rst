Running the test set
====================


How to run tests using runtest
------------------------------

It is *very important* that you verify that your DIRAC installation correctly
reproduces the reference test set *before* running any production calculations.

To get documentation on the test set, run::

  ./runtest --help

With the following command you run the entire test set::

  ./runtest --all

If everything goes well the script will report calculations completed
and provide a summary (number of tests may be different)::

  summary:   114 ok

  total time: 00h03m42s

If there is a problem the script will report that one or more tests failed
(wrong output) or crashed (execution error).  You may also run individual tests
with::

  ./runtest --tests="test1 test2"

where test1 and test2 should be replaced by the name(s) of the test(s)
that you want to rerun.

If you are in the subdirectory of an individual test (example: dft_ac), you can
do this instead to run just this one test::

  $ [cd test/dft_ac]
  $ ../../runtest


How to run tests using CTest
----------------------------

Go to the build directory::

  $ cd build

Now run all tests::

  $ make test

You can also run individual tests (example: dft_ac)::

  $ ctest -R dft_ac

The "-R" means regular expression! So you can run all tests
that start with "dft" using::

  $ ctest -R dft

By the way behind the curtains CTest uses runtest.


Some tests are broken, what now?
--------------------------------

You are a user? Then please inform the developers via the mailing list.

You are a developer? Perhaps it is not your fault. Have a look at the nightly
testing dashboard where you can check whether this test is marked as red
(broken) or green (passing). If it is green there but broken on your machine
it can either mean that the test has trouble with your architecture/compiler
(in this case inform other developers) or that you broke it with recent changes
- then you have to go in and fix it before you git push to origin/master.


Suddenly some red tests appear on the nightly testing dashboard, what happens now?
----------------------------------------------------------------------------------

Miro or Radovan or somebody else examining the nightly testing dashboard will
notice it and write to the developer who is likely to have introduced the
conflicting change.


