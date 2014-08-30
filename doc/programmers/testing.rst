

Never commit functionality to the main development line without tests
---------------------------------------------------------------------

If you commit functionality to the main development line without tests then
this functionality will break sooner or later and we have no automatic
mechanism to detect it. Committing new code without tests is bad karma.


When you add a new test always add it to pamadm
-----------------------------------------------

Otherwise the test will never be run as part of the default test suite.


Strive for portability
----------------------

Avoid shell programming or symlinks in test scripts otherwise the tests are not
portable to Windows. Therefore do not use os.system() or os.symlink(). Do not
use explicit forward slashes for paths, instead use os.path.join().


Never add inputs to the test directories which are never run
------------------------------------------------------------

We want all inputs and outputs inside test/ to be accessile by the default test
suite. Otherwise we have no automatic way to detect that some inputs or outputs
have degraded. And degraded inputs and outputs are useless and confusing.


How can I run a single test?
----------------------------

The test runner is CTest. You can run a single test with::

  ctest -R mynewtest

Alternatively you can run a test directly in the test directory.
For this check out::

  cd test/mynewtest
  ./test --help


Always test that the test really works
--------------------------------------

It is easy to make a mistake and create a test which is always "successful".
Test that your test catches mistakes. Verify whether it extracts the right
numbers.


How to add labels to tests
--------------------------

Here is an example::

  add_test('dft_ac', 'short;dft')

If you give several labels, separate them with ";".


How to add a new test
---------------------

All you need to add a test is to add a python script called "test" in a new
directory.

In most cases you probably want to use the runtest_dirac.py test library
since it gives you the framework to run DIRAC jobs and
to filter numbers with numerical tolerance.  But you don't have to. You have
full liberty to test whatever you want in whatever way you want. Important is
that you return an exit code. If the exit code is 0, the test was successful,
if non-0, then the test has failed. Add the test to pamadm also.
