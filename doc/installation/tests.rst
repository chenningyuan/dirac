

Running the test set
====================

It is *very important* that you verify that your DIRAC installation correctly
reproduces the reference test set *before* running any production calculations.

To get documentation on the test set, run::

  ./runtest --help

As a starter try the quick test that covers most of the standard run
types::

  ./runtest --quick

If everything goes well the script will report calculations completed
and provide a summary::

  summary:   114 ok

  total time: 00h03m42s

If there is a problem the script will report that one or more tests failed
(wrong output) or crashed (execution error).  You may also run individual tests
with::

  ./runtest --tests="test1 test2"

where test1 and test2 should be replaced by the name(s) of the test(s)
that you want to rerun. When everything seems fine try also the complete
test set::

  ./runtest --all
